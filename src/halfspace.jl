"""
    StaticNgon(point_type[, capacity=32])
    StaticNgon(polygon[, capacity])

Fixed-capacity mutable polygon used by the allocation-free `intersect!`,
`shift`, and `reconstruct!` paths. `capacity` is the number of vertices that
may be written; choose it large enough for every intermediate polygon.
"""
mutable struct StaticNgon{N, P}
    const vertices::MVector{N, P}
    nr_verts::Int
    interface_index::Int
end

StaticNgon(P, N::Integer=32) = StaticNgon{N, P}(MVector{N, P}(undef), 0, 0)
StaticNgon(p::Ngon, N::Integer=max(8, 2length(p.vertices) + 2)) =
    StaticNgon(eltype(p.vertices), N)

"""Maximum number of vertices that a `StaticNgon` can store."""
capacity(poly::StaticNgon) = length(poly.vertices)

abstract type HalfSpace{D} end

import Meshes: Ngon

"""
    PlanarHS(normal, shift)
    PlanarHS(angle, area, polygon; workspace, shift_workspace)

A two-dimensional half-space satisfying `normal ⋅ x ≤ shift`. Normals are
stored as immutable `SVector`s, including when constructed from ordinary
vectors. The three-argument form computes the shift that encloses `area` of
`polygon`; `normal` is generated from `angle` and is therefore unit length.
"""
struct PlanarHS{D, T<:Real, Q<:Quantity} <: HalfSpace{D}
    𝛈::SVector{D, T}
    shift::Q
end

function PlanarHS(𝛈::AbstractVector{T}, shift::Q) where {T<:Real, Q<:Quantity}
    D = length(𝛈)
    return PlanarHS{D, T, Q}(SVector{D, T}(𝛈), shift)
end

function PlanarHS{D}(𝛈::AbstractVector{T}, shift::Q) where {D, T<:Real, Q<:Quantity}
    return PlanarHS{D, T, Q}(SVector{D, T}(𝛈), shift)
end

"""Return the outward normal of a planar half-space."""
normal(p::PlanarHS) = p.𝛈

complement(p::PlanarHS) = PlanarHS(-p.𝛈, -p.shift)

@inline distance(p::PlanarHS{2}, x::Point) =
    p.𝛈[1] * x.coords.x + p.𝛈[2] * x.coords.y - p.shift

function PlanarHS(θ::T, αvol::Quantity, c::Ngon;
    workspace::StaticNgon=StaticNgon(c),
    shift_workspace::AbstractVector{<:Real}=Vector{Float64}(undef, length(c.vertices)),
) where {T<:Real}
    𝛈 = angle_to_normal(θ)
    return PlanarHS{2}(𝛈, shift(c, 𝛈, αvol;
        workspace=workspace, shift_workspace=shift_workspace))
end

function PlanarHS(θ::T, αvol::Quantity, c::StaticNgon{N, P};
    workspace::StaticNgon=StaticNgon(P, N + 2),
    shift_workspace::AbstractVector{<:Real}=MVector{N + 2, Float64}(undef),
) where {T<:Real, N, P<:Point}
    𝛈 = angle_to_normal(θ)
    return PlanarHS{2}(𝛈, shift(c, 𝛈, αvol;
        workspace=workspace, shift_workspace=shift_workspace))
end

function PlanarHS(v1::Point, v2::Point)
    𝛈 = angle_to_normal(atan(v1.coords.x - v2.coords.x, v2.coords.y - v1.coords.y))
    return PlanarHS{2}(𝛈, 𝛈 ⋅ to(v1))
end

import Base: intersect

"""
    intersect(polygon, halfspace)

Return the part of a 2D polygon inside `halfspace`, or `nothing` when it is
empty. Use `intersect!` with a `StaticNgon` when this is a hot path.
"""
Base.intersect(c::Ngon, p::PlanarHS{2}; kwargs...) =
    Ngon(Base.intersect!(StaticNgon(c), c, p; kwargs...))

Ngon(poly::StaticNgon) = poly.nr_verts < 3 ? nothing : Ngon(poly.vertices[1:poly.nr_verts]...)

function _require_capacity(out::StaticNgon, required::Integer)
    capacity(out) ≥ required || throw(ArgumentError(
        "StaticNgon capacity $(capacity(out)) is too small; clipping needs at least $required vertices"))
    return nothing
end

function Base.intersect!(out::StaticNgon{N, P}, c::Ngon, p::PlanarHS{2};
    tol::Real=√eps(typeof(c.vertices[1].coords.x.val)),
) where {N, P<:Point}
    intersect!(out, c.vertices, p; tol=tol)
end

function Base.intersect!(out::StaticNgon{N, P}, input::StaticNgon{M, P}, p::PlanarHS{2};
    tol::Real=√eps(),
) where {N, M, P<:Point}
    intersect!(out, view(input.vertices, 1:input.nr_verts), p; tol=tol)
end

function Base.intersect!(out::StaticNgon{N, P}, verts::AbstractVector{P}, p::PlanarHS{2};
    tol::Real=√eps(),
) where {N, P<:Point}
    nr_old_verts = length(verts)
    nr_old_verts == 0 && return (out.nr_verts = 0; out.interface_index = 0; out)
    _require_capacity(out, nr_old_verts + 1)

    out.nr_verts = 0
    out.interface_index = 0
    any_bisected = false

    next_dist = zero(distance(p, verts[1]))
    next_inside = false
    for (index, curr_vert) in enumerate(verts)
        next_index = mod1(index + 1, nr_old_verts)
        if index == 1
            curr_dist = distance(p, curr_vert)
            curr_inside = curr_dist ≤ zero(curr_dist)
        else
            curr_dist = next_dist
            curr_inside = next_inside
        end
        next_vert = verts[next_index]
        next_dist = distance(p, next_vert)
        next_inside = next_dist ≤ zero(next_dist)

        if curr_inside
            out.nr_verts += 1
            out.vertices[out.nr_verts] = curr_vert
        end

        if curr_inside != next_inside
            any_bisected = true
            coeff = abs(curr_dist / (next_dist - curr_dist))
            if (curr_inside && coeff > tol) || (next_inside && coeff < 1 - tol)
                out.nr_verts += 1
                out.vertices[out.nr_verts] = curr_vert + coeff * (next_vert - curr_vert)
            end
            if out.interface_index == 0
                out.interface_index = curr_inside ? out.nr_verts : out.nr_verts - 1
            end
        end
    end

    if any_bisected && out.nr_verts > 0
        out.interface_index = mod1(out.interface_index, out.nr_verts)
    end
    return out
end

function copy!(to::StaticNgon{N1, P}, from::StaticNgon{N2, P}) where {N1, N2, P<:Point}
    _require_capacity(to, from.nr_verts)
    to.nr_verts = from.nr_verts
    to.interface_index = from.interface_index
    for i in 1:from.nr_verts
        to.vertices[i] = from.vertices[i]
    end
    return to
end

function copy!(to::StaticNgon{N, P}, from::Ngon) where {N, P<:Point}
    _require_capacity(to, length(from.vertices))
    to.nr_verts = length(from.vertices)
    to.interface_index = 0
    for i in eachindex(from.vertices)
        to.vertices[i] = from.vertices[i]
    end
    return to
end

function Base.intersect!(out::StaticNgon{N, P}, c1::Ngon, c2::Ngon;
    tol::Real=√eps(typeof(c1.vertices[1].coords.x.val)),
    workspace::StaticNgon=StaticNgon(c1),
) where {N, P<:Point}
    static_c1 = StaticNgon(c1)
    copy!(static_c1, c1)
    return intersect!(out, static_c1, c2; tol=tol, workspace=workspace)
end

function Base.intersect!(out::StaticNgon{N1, P}, c1::StaticNgon{N2, P}, c2::Ngon;
    tol::Real=√eps(typeof(c1.vertices[1].coords.x.val)),
    workspace::StaticNgon=StaticNgon(c2),
) where {N1, N2, P<:Point}
    copy!(workspace, c1)
    for (index, curr_vert) in enumerate(c2.vertices)
        next_vert = c2.vertices[mod1(index + 1, length(c2.vertices))]
        intersect!(out, workspace, PlanarHS(curr_vert, next_vert); tol=tol)
        index < length(c2.vertices) && copy!(workspace, out)
    end
    return out
end

"""Signed area of a fixed-capacity polygon."""
function smeasure(poly::StaticNgon{N, P}) where {N, P<:Point}
    poly.nr_verts == 0 && return 0u"m^2"
    poly.nr_verts < 3 && return zero(poly.vertices[1].coords.x * poly.vertices[1].coords.y)
    area = zero(poly.vertices[1].coords.x * poly.vertices[1].coords.y)
    previous = poly.vertices[poly.nr_verts]
    for current in view(poly.vertices, 1:poly.nr_verts)
        area += previous.coords.x * current.coords.y - previous.coords.y * current.coords.x
        previous = current
    end
    return area / 2
end

"""
    smeasure(halfspace, polygon; workspace)

Area of the portion of `polygon` inside `halfspace`. Supply a reusable
`StaticNgon` workspace to avoid allocations in a loop.
"""
function smeasure(p::PlanarHS{2}, c::Ngon; workspace::StaticNgon=StaticNgon(c))
    intersect!(workspace, c, p)
    return workspace.nr_verts < 3 ?
        zero(c.vertices[1].coords.x * c.vertices[1].coords.y) : smeasure(workspace)
end

function smeasure(p::PlanarHS{2}, c::StaticNgon;
    workspace::StaticNgon=StaticNgon(eltype(c.vertices), capacity(c) + 2),
)
    intersect!(workspace, c, p)
    return workspace.nr_verts < 3 ? smeasure(c) * 0 : smeasure(workspace)
end

"""
    shift(polygon, normal, area; workspace, shift_workspace)

Find the half-space shift whose clipped signed area is `area`. `area` must lie
between zero and the polygon's signed area. The default is safe for polygons
of any size; pass reusable workspaces in repeated calls.
"""
function shift(c::Ngon, 𝛈::AbstractVector, αvol::Quantity;
    workspace::StaticNgon=StaticNgon(c),
    shift_workspace::AbstractVector{<:Real}=Vector{Float64}(undef, length(c.vertices)),
)
    return _shift(c, SVector{2}(𝛈), αvol, workspace, shift_workspace)
end

function shift(c::StaticNgon{N, P}, 𝛈::AbstractVector, αvol::Quantity;
    workspace::StaticNgon=StaticNgon(P, N + 2),
    shift_workspace::AbstractVector{<:Real}=MVector{N + 2, Float64}(undef),
) where {N, P<:Point}
    return _shift(c, SVector{2}(𝛈), αvol, workspace, shift_workspace)
end

function _shift(c, 𝛈::SVector{2}, αvol::Quantity, workspace::StaticNgon,
    shift_workspace::AbstractVector{<:Real})
    nverts = c isa Ngon ? length(c.vertices) : c.nr_verts
    nverts > 0 || throw(ArgumentError("cannot find a shift for an empty polygon"))
    length(shift_workspace) ≥ nverts || throw(ArgumentError(
        "shift_workspace has length $(length(shift_workspace)); at least $nverts values are required"))
    _require_capacity(workspace, nverts + 1)

    c_measure = smeasure(c)
    zero_area = zero(c_measure)
    zero_area ≤ αvol ≤ c_measure || throw(DomainError(αvol,
        "area must lie between zero and the polygon's signed area"))

    verts = c isa Ngon ? c.vertices : view(c.vertices, 1:c.nr_verts)
    shift_unit = oneunit(𝛈[1] * verts[1].coords.x + 𝛈[2] * verts[1].coords.y)
    for index in eachindex(verts)
        vertex = verts[index]
        shift_workspace[index] = ustrip((𝛈[1] * vertex.coords.x + 𝛈[2] * vertex.coords.y) / shift_unit)
    end
    sort!(view(shift_workspace, 1:nverts))

    αvol == zero_area && return shift_workspace[1] * shift_unit
    αvol == c_measure && return shift_workspace[nverts] * shift_unit

    αerr(shift_value) = smeasure(intersect!(workspace, c,
        PlanarHS{2}(𝛈, shift_value * shift_unit))) - αvol

    αerr_previous = -αvol
    for index in 2:nverts
        αerr_current = index == nverts ? c_measure - αvol : αerr(shift_workspace[index])
        αerr_current == zero_area && return shift_workspace[index] * shift_unit

        if sign(αerr_current) != sign(αerr_previous)
            shift0 = shift_workspace[index - 1]
            shift2 = shift_workspace[index]
            shift1 = (shift0 + shift2) / 2
            h = shift2 - shift1
            h == 0 && continue

            αerr0 = ustrip(αerr_previous / (shift_unit * shift_unit))
            αerr1 = ustrip(αerr(shift1) / (shift_unit * shift_unit))
            αerr2 = ustrip(αerr_current / (shift_unit * shift_unit))
            A = (0.5αerr2 - αerr1 + 0.5αerr0) / h^2
            B = (αerr2 - αerr0) / (2h)
            C = αerr1
            x1, x2, _ = parabola_roots(A, B, C)
            root = shift0 ≤ x1 + shift1 ≤ shift2 ? x1 : x2
            return (shift1 + root) * shift_unit
        end
        αerr_previous = αerr_current
    end

    error("GeometricVOF.shift: no bracket found — possible numerical degeneracy (area=$αvol)")
end

function shift_extrema(c::Ngon, 𝛈::AbstractVector)
    verts = c.vertices
    first_shift = 𝛈[1] * verts[1].coords.x + 𝛈[2] * verts[1].coords.y
    shift_min = first_shift
    shift_max = first_shift
    for vertex in @view verts[2:end]
        value = 𝛈[1] * vertex.coords.x + 𝛈[2] * vertex.coords.y
        value < shift_min && (shift_min = value)
        value > shift_max && (shift_max = value)
    end
    return shift_min, shift_max
end

"""
    sorted_unique_approx(polygon; tol)

Remove adjacent approximately-equal vertices, returning `nothing` if fewer
than three remain.
"""
function sorted_unique_approx(c::Ngon; tol::Real=√eps(typeof(c.vertices[1].coords.x.val)))
    vs = vertices(c)
    remove = Int[]
    for (index, v1) in enumerate(vs)
        v2 = vs[index == length(vs) ? 1 : index + 1]
        if ustrip(abs(v1.coords.x - v2.coords.x)) < tol &&
           ustrip(abs(v1.coords.y - v2.coords.y)) < tol
            push!(remove, index)
        end
    end
    isempty(remove) && return c
    length(vs) - length(remove) < 3 && return nothing
    return Ngon(vs[setdiff(eachindex(vs), remove)]...)
end
