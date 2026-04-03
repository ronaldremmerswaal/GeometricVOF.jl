# const PlanarHS = Tuple{Vector{T}, Quantity} where {T<:Real}

mutable struct StaticNgon{N, P}
    const vertices::MVector{N, P}
    nr_verts::Int64
    interface_index::Int64
end
StaticNgon(P, N::Int=32) = StaticNgon{N, P}(MVector{N, P}(undef), 0, 0)
StaticNgon(p::Ngon, N::Int=32) = StaticNgon(eltype(p.vertices), N)

abstract type HalfSpace{D} end

"""
    PlanarHS(𝛈, shift)

Represents the half-space defined by

    𝛈 ⋅ 𝐱 ≤ shift
"""
struct PlanarHS{D, V, Q} <: HalfSpace{D}
    𝛈::V
    shift::Q
end
# PlanarHS(𝛈::Vector, shift::Number) = PlanarHS{length(𝛈)}(SVector{length(𝛈)}(𝛈), shift * u"m")
PlanarHS(𝛈::V, shift::Q) where {V <: AbstractVector, Q <: Quantity} = PlanarHS{length(𝛈), V, Q}(SVector{length(𝛈)}(𝛈), shift)
PlanarHS{D}(𝛈::V, shift::Q) where {D, V <: AbstractVector, Q <: Quantity} = PlanarHS{D, V, Q}(SVector{D}(𝛈), shift)

complement(p::PlanarHS) = PlanarHS(-p.𝛈, -p.shift)

distance(p::PlanarHS, 𝐱::Point) = p.𝛈[1] * 𝐱.coords.x + p.𝛈[2] * 𝐱.coords.y - p.shift

function PlanarHS(θ::T, αvol::Quantity, c::Ngon; workspace::StaticNgon=StaticNgon(c), shift_workspace::MVector=MVector{32, Float64}(undef)) where {T <: Real}
    𝛈 = GeometricVOF.angle_to_normal(θ)
    s = shift(c, 𝛈, αvol; workspace=workspace, shift_workspace=shift_workspace)
    return PlanarHS{2}(𝛈, s)
end

function PlanarHS(v1::Point, v2::Point)
    𝛈 = GeometricVOF.angle_to_normal(atan(v1.coords.x - v2.coords.x, v2.coords.y - v1.coords.y))
    return PlanarHS{2}(𝛈, 𝛈 ⋅ to(v1))
end

import Base: intersect
"""
    intersect(c, p)

Intersects the polygon `c` with halfspace defined by `p`.

# Examples
```julia-repl
julia> c = Triangle((0., 0.), (1., 0.), (0., 1.))
julia> p = PlanarHS([-1., 0.], -0.5)
julia> intersect(c, p)
Triangle
├─ Point(x: 0.5 m, y: 0.5 m)
├─ Point(x: 0.5 m, y: 0.0 m)
└─ Point(x: 1.0 m, y: 0.0 m)
```
"""

Base.intersect(c::Ngon, p::PlanarHS{2}; kwargs...) =
    Ngon(Base.intersect!(StaticNgon(eltype(c.vertices)), c, p; kwargs...))

Ngon(poly::StaticNgon) = poly.nr_verts < 3 ? nothing : Ngon(poly.vertices[1:poly.nr_verts]...)
function smeasure(p::PlanarHS, c::Ngon)
    cp = intersect(c, p)
    if isnothing(cp)
        return 0u"m^2"
    else
        return smeasure(cp)
    end
end

function Base.intersect!(out::StaticNgon{N, P}, c::Ngon, p::PlanarHS{2}; tol::Real=√eps(typeof(c.vertices[1].coords.x.val))) where {N, P<:Point}
    intersect!(out, c.vertices, p; tol=tol)
end

function Base.intersect!(out::StaticNgon{N, P}, in::StaticNgon{N, P}, p::PlanarHS{2}; tol::Real=√eps(typeof(c.vertices[1].coords.x.val))) where {N, P<:Point}
    intersect!(out, view(in.vertices, 1:in.nr_verts), p; tol=tol)
end

function Base.intersect!(out::StaticNgon{N, P}, verts::AbstractVector{P}, p::PlanarHS{2}; tol::Real=√eps()) where {N, P<:Point}

    nr_old_verts = length(verts)

    # Construct new polygon by looping over the edges of the old polygon
    out.nr_verts = 0
    next_dist = 0u"m"
    next_inside = false
    out.interface_index = 0
    any_bisected = false
    for (cdx, curr_vert) ∈ enumerate(verts)
        ndx = mod1(cdx + 1, nr_old_verts)

        if cdx == 1
            curr_dist = distance(p, curr_vert)
            curr_inside = curr_dist ≤ 0u"m"
        else
            curr_dist = next_dist
            curr_inside = next_inside
        end
        next_vert = verts[ndx]
        next_dist = distance(p, next_vert)
        next_inside = next_dist ≤ 0u"m"

        if curr_inside
            out.nr_verts += 1
            out.vertices[out.nr_verts] = curr_vert
        end

        edge_is_bisected = curr_inside != next_inside
        any_bisected = any_bisected || edge_is_bisected

        if edge_is_bisected
            coeff = abs(curr_dist / (next_dist - curr_dist))
            if (curr_inside && coeff > tol) ||
               (next_inside && coeff < 1 - tol)
                out.nr_verts += 1
                out.vertices[out.nr_verts] = curr_vert + coeff * (next_vert - curr_vert)
            end
            if out.interface_index == 0
                if curr_inside
                    out.interface_index = out.nr_verts
                else
                    out.interface_index = out.nr_verts - 1
                end
            end
        end
    end

    if any_bisected
        out.interface_index = mod1(out.interface_index, out.nr_verts)
    end

    return out
end

function copy!(to::StaticNgon{N1, P}, from::StaticNgon{N2, P}) where {N1, N2, P<:Point}
    to.nr_verts = from.nr_verts
    to.interface_index = from.interface_index
    for i ∈ 1:from.nr_verts
        to.vertices[i] = from.vertices[i]
    end
    return to
end

function copy!(to::StaticNgon{N, P}, from::Ngon) where {N, P<:Point}
    to.nr_verts = length(from.vertices)
    to.interface_index = 0
    for i ∈ 1:to.nr_verts
        to.vertices[i] = from.vertices[i]
    end
    return to
end

function Base.intersect!(out::StaticNgon{N, P}, c1::Ngon, c2::Ngon;
    tol::Real=√eps(typeof(c1.vertices[1].coords.x.val)), workspace::StaticNgon=StaticNgon(c1))  where {N, P<:Point}

    static_c1 = StaticNgon(c1)
    copy!(static_c1, c1)

    intersect!(out, static_c1, c2; tol=tol, workspace=workspace)
end

function Base.intersect!(out::StaticNgon{N1, P}, c1::StaticNgon{N2, P}, c2::Ngon;
    tol::Real=√eps(typeof(c1.vertices[1].coords.x.val)), workspace::StaticNgon=StaticNgon(c1))  where {N1, N2, P<:Point}

    copy!(workspace, c1)

    for (cdx, curr_vert) ∈ enumerate(c2.vertices)
        next_vert = c2.vertices[mod1(cdx + 1, length(c2.vertices))]
        hs = PlanarHS(curr_vert, next_vert)

        intersect!(out, workspace, hs; tol=tol)

        if cdx < length(c2.vertices)
            copy!(workspace, out)
        end
    end

    return out
end

function smeasure(poly::StaticNgon{N, P}) where {N, P<:Point}
    M = 0u"m^2"
    if poly.nr_verts < 3
        return M
    end

    for vdx = 1:poly.nr_verts
        v1 = poly.vertices[vdx]
        v2 = poly.vertices[mod1(vdx + 1, poly.nr_verts)]
        M += v1.coords.x * v2.coords.y - v1.coords.y * v2.coords.x
    end
    M /= 2
end


"""
    shift(c, 𝛈, α)

Shift such that shifted plane with normal `𝛈` yields intersection volume given by
`α`. It is the inverse of the `measure(c, Plane(𝛈, shift))` function.

# Examples
```julia-repl
julia> c = Triangle((0., 0.), (1., 0.), (0., 1.))
julia> shift(c, [1.0, 0.0], 0.21875u"m^2")
0.25u"m"
```
"""
shift(c::Ngon, 𝛈::Vector, α::Quantity) = shift(c, SVector{2}(𝛈), α)
function shift(c::Ngon, 𝛈::SVector{2}, αvol::Quantity; workspace::StaticNgon=StaticNgon(c), shift_workspace::MVector=MVector{32, Float64}(undef))

    α_err(p) = smeasure(intersect!(workspace, c, p)) - αvol

    # Determine shift at each vertex
    n_verts = length(c.vertices)
    for (vdx, v) ∈ enumerate(c.vertices)
        shift_workspace[vdx] = ustrip(𝛈[1] * v.coords.x + 𝛈[2] * v.coords.y)
    end

    # Evaluate the error function at the shift_workspace
    sort!(view(shift_workspace, 1:n_verts))
    α_err_prev = -αvol

    c_measure = smeasure(c)

    if αvol == 0u"m^2"
        return shift_workspace[1]u"1m"
    elseif αvol == c_measure
        return shift_workspace[n_verts]u"1m"
    end

    for i ∈ 2:n_verts
        if i == n_verts
            α_err_curr = c_measure - αvol
        else
            α_err_curr = α_err(PlanarHS{2}(𝛈, shift_workspace[i]u"1m"))
        end

        if α_err_curr == 0u"m^2" return shift_workspace[i]u"1m" end

        if sign(α_err_curr) != sign(α_err_prev)
            # We have found a sign change, so we can use the two shift_workspace to bracket the root
            shift0 = shift_workspace[i - 1]u"1m"
            shift2 = shift_workspace[i]u"1m"
            shift1 = (shift0 + shift2) / 2

            # Moreover, the dependence in the bracket is quadratic, so 3 values are sufficient
            α_err0 = ustrip(α_err_prev)
            α_err1 = ustrip(α_err(PlanarHS{2}(𝛈, shift1)))
            α_err2 = ustrip(α_err_curr)

            h = ustrip(shift2 - shift1)
            A = (.5α_err2 - α_err1 + .5α_err0) / h^2
            B = (α_err2 - α_err0) / 2h
            C = α_err1

            x1, x2, = parabola_roots(A, B, C)
            x1u = shift1 + x1*u"m"
            x2u = shift1 + x2*u"m"

            return shift0 ≤ x1u ≤ shift2 ? x1u : x2u
        end

        α_err_prev = α_err_curr
    end

    # Unreachable for valid inputs (monotone bracket guaranteed by IVT).
    # Explicit error keeps the return type as Quantity (not Union{Quantity,Nothing}).
    error("GeometricVOF.shift: no bracket found — possible numerical degeneracy (αvol=$αvol)")
end

function shift_extrema(c::Ngon, 𝛈::SVector{2})
    shift_min = floatmax()u"m"
    shift_max = floatmin()u"m"
    for v ∈ c.vertices
        shift_val = 𝛈[1] * v.coords.x + 𝛈[2] * v.coords.y
        if shift_val < shift_min
            shift_min = shift_val
        end
        if shift_val > shift_max
            shift_max = shift_val
        end
    end
    return shift_min, shift_max
end

"""
    sorted_unique_approx(c::Ngon; tol::Real)

Remove subsequent vertices that are approximately equal to each other.
"""
function sorted_unique_approx(c::Ngon; tol::Real=√eps(typeof(c.vertices[1].coords.x.val)))
    vs = vertices(c)
    rm_indices = Vector{Int}()

    for (vdx, v1) ∈ enumerate(vs)
        v2 = vs[vdx == length(vs) ? 1 : vdx + 1]

        if ustrip(abs(v1.coords.x - v2.coords.x)) < tol && ustrip(abs(v1.coords.y - v2.coords.y)) < tol
            push!(rm_indices, vdx)
        end
    end

    if isempty(rm_indices)
        return c
    elseif length(vs) - length(rm_indices) < 3
        return nothing
    else
        return Ngon(vs[setdiff(1:length(vs), rm_indices)]...)
    end

end
