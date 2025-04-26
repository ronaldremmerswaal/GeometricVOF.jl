abstract type HalfSpace{D} end

"""
    PlanarHS(𝛈, shift)

Represents the half-space defined by

    𝛈 ⋅ 𝐱 ≤ shift
"""
struct PlanarHS{D} <: HalfSpace{D}
    𝛈::SVector{D}
    shift::Quantity
end
# PlanarHS(𝛈::Vector, shift::Number) = PlanarHS{length(𝛈)}(SVector{length(𝛈)}(𝛈), shift * u"m")
PlanarHS(𝛈::Vector, shift::Quantity) = PlanarHS{length(𝛈)}(SVector{length(𝛈)}(𝛈), shift)

complement(p::PlanarHS) = PlanarHS(-p.𝛈, -p.shift)

distance(p::PlanarHS{2}, 𝐱::Point) = p.𝛈[1] * 𝐱.coords.x + p.𝛈[2] * 𝐱.coords.y - p.shift

function PlanarHS(θ::Real, αvol::Quantity, c::Ngon)
    𝛈 = GeometricVOF.angle_to_normal(θ)
    s = shift(c, 𝛈, αvol)
    return PlanarHS{2}(𝛈, s)
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
function Base.intersect(c::Ngon, p::PlanarHS{2})

    nr_old_verts = length(vertices(c))

    # Compute the distance of each vertex to the plane
    dist = zeros(nr_old_verts)u"m"
    nr_inside = 0
    for (vdx, vert) ∈ enumerate(c.vertices)
        dist[vdx] = distance(p, vert)
        if dist[vdx] ≤ 0u"m"
            nr_inside += 1
        end
    end

    # Quick return in trivial cases
    if nr_inside == 0
        return nothing  # TODO type stability
    elseif nr_inside == length(c.vertices)
        return c        # TODO no copy available
    end

    # Find first bisected edge; this ensures that the first edge of the new polygon coincides
    # with the HalfSpace
    start_index = 0
    for (vdx, d) ∈ enumerate(dist)
        ndx = vdx == nr_old_verts ? 1 : vdx + 1
        nd = dist[ndx]

        edge_is_bisected = (d ≤ 0u"m") != (nd ≤ 0u"m")
        if edge_is_bisected
            start_index = vdx
            break
        end
    end

    new_verts = Vector{eltype(c.vertices)}()

    # Construct new polygon by looping over the edges of the old polygon
    vdx = start_index
    for _ ∈ 1:nr_old_verts
        ndx = vdx == nr_old_verts ? 1 : vdx + 1
        edge_is_bisected = (dist[vdx] ≤ 0u"m") != (dist[ndx] ≤ 0u"m")

        if edge_is_bisected
            coeff = abs(dist[vdx] / (dist[ndx] - dist[vdx]))
            if (coeff < 1 || dist[ndx] > 0u"m") && (coeff > 0 || dist[vdx] > 0u"m")
                push!(new_verts, c.vertices[vdx] +
                    coeff * (c.vertices[ndx] - c.vertices[vdx]))
            end
        end

        if dist[ndx] ≤ 0u"m"
            push!(new_verts, c.vertices[ndx])
        end

        vdx = ndx
    end

    if length(new_verts) < 3
        return nothing
    else
        return Ngon(new_verts...)
    end
end

function measure(p::PlanarHS{2}, c::Ngon)
    cp = c ∩ p
    isnothing(cp) ? 0u"m^2" : measure(cp)
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
function shift(c::Ngon, 𝛈::SVector{2}, α::Quantity) # TODO constrain α to have units m^2

    α_err(shift) = measure(PlanarHS(𝛈, shift), c) - α

    # Determine shift at each vertex
    shifts = Vector{Quantity}(undef, length(c.vertices))
    for (vdx, v) ∈ enumerate(c.vertices)
        shifts[vdx] =  𝛈[1] * v.coords.x + 𝛈[2] * v.coords.y
    end

    # Evaluate the error function at the shifts
    perm = sortperm(shifts)
    shifts = shifts[perm]
    α_err_prev = -α
    if α_err_prev == 0u"m^2" return shifts[1] end

    for i ∈ 2:length(c.vertices)
        if i < length(c.vertices)
            α_err_curr = α_err(shifts[i])
        else
            α_err_curr = measure(c) - α
        end

        if α_err_curr == 0u"m^2" return shifts[i] end

        if sign(α_err_curr) != sign(α_err_prev)
            # We have found a sign change, so we can use the two shifts to bracket the root
            shift0 = shifts[i - 1]
            shift2 = shifts[i]
            shift1 = (shift0 + shift2) / 2

            # Moreover, the dependence in the bracket is quadratic, so 3 values are sufficient
            α_err0 = ustrip(α_err_prev)
            α_err1 = ustrip(α_err(shift1))
            α_err2 = ustrip(α_err_curr)

            h = ustrip(shift2 - shift1)
            A = (.5α_err2 - α_err1 + .5α_err0) / h^2
            B = (α_err2 - α_err0) / 2h
            C = α_err1

            rts = shift1 .+ roots(Polynomial([C, B, A]))u"m"

            return shift0 ≤ rts[1] ≤ shift2 ? rts[1] : rts[2]
        end

        α_err_prev = α_err_curr
    end

    # bracket = α_err0 > 0u"m^2" ? (shift_min, shift0) : (shift0, shift_max)


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
