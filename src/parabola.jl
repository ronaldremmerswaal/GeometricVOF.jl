"""
    Parabola(normal, shift, curvature, origin)

An oriented parabolic interface in two dimensions.  Its liquid side is
``normal ⋅ (x - origin) - shift + curvature / 2 * (tangent ⋅ (x - origin))^2 ≤ 0``.
`normal` is unit length, `curvature` has inverse-length units, and `origin`
is the point about which the local parabolic coordinates are measured.
"""
struct Parabola{T<:Real,S<:Quantity,K<:Quantity,P<:Point}
    𝛈::SVector{2,T}
    shift::S
    curvature::K
    origin::P
end

function Parabola(𝛈::AbstractVector{T}, shift::S, curvature::K, origin::P) where
    {T<:Real,S<:Quantity,K<:Quantity,P<:Point}
    length(𝛈) == 2 || throw(DimensionMismatch("a 2D parabola needs a two-component normal"))
    return Parabola{T,S,K,P}(SVector{2,T}(𝛈), shift, curvature, origin)
end

normal(p::Parabola) = p.𝛈
tangent(p::Parabola) = SVector(-p.𝛈[2], p.𝛈[1])
complement(p::Parabola) = Parabola(-p.𝛈, -p.shift, -p.curvature, p.origin)

@inline function distance(p::Parabola, x::Point)
    dx = x.coords.x - p.origin.coords.x
    dy = x.coords.y - p.origin.coords.y
    τ = -p.𝛈[2] * dx + p.𝛈[1] * dy
    return p.𝛈[1] * dx + p.𝛈[2] * dy - p.shift + p.curvature * τ^2 / 2
end

"""
    StaticParabolicNgon(polygon, parabola[, capacity])

Fixed-capacity polygonal representation of a clipped parabolic region.  Each
entry in `parabolic_faces` describes the edge beginning at the corresponding
vertex.  A `true` value means that edge is the parabolic arc, rather than the
straight chord joining its endpoints.  This keeps clipping and exact moment
evaluation allocation-free after workspace construction.
"""
mutable struct StaticParabolicNgon{N,P,C,M}
    const vertices::MVector{N,P}
    const parabolic_faces::MVector{N,Bool}
    nr_verts::Int
    curve::C
    complement::Bool
    original_moments::M
end

function StaticParabolicNgon(P, p::Parabola, original_moments, N::Integer=32)
    return StaticParabolicNgon{N,P,typeof(p),typeof(original_moments)}(
        MVector{N,P}(undef), MVector{N,Bool}(undef), 0, p, false, original_moments)
end

function StaticParabolicNgon(c::Ngon, p::Parabola,
    N::Integer=max(8, 2length(c.vertices) + 2))
    return StaticParabolicNgon(eltype(c.vertices), p, moments(c), N)
end

function StaticParabolicNgon(c::StaticNgon{N,P}, p::Parabola,
    capacity::Integer=max(8, 2c.nr_verts + 2)) where {N,P}
    return StaticParabolicNgon(P, p, moments(c), capacity)
end

"""Maximum number of vertices that a `StaticParabolicNgon` can store."""
capacity(poly::StaticParabolicNgon) = length(poly.vertices)

function _require_capacity(out::StaticParabolicNgon, required::Integer)
    capacity(out) ≥ required || throw(ArgumentError(
        "StaticParabolicNgon capacity $(capacity(out)) is too small; clipping needs at least $required vertices"))
    return nothing
end

@inline function _polygon_moments(vertices, nverts::Integer, fallback)
    nverts < 3 && return zero(fallback[1]), zero(fallback[2])
    previous = vertices[nverts]
    area = zero(previous.coords.x * previous.coords.y)
    moment_x = zero(previous.coords.x * previous.coords.x * previous.coords.y)
    moment_y = zero(previous.coords.x * previous.coords.y * previous.coords.y)
    for index in 1:nverts
        current = vertices[index]
        cross = previous.coords.x * current.coords.y - previous.coords.y * current.coords.x
        area += cross
        moment_x += cross * (previous.coords.x + current.coords.x)
        moment_y += cross * (previous.coords.y + current.coords.y)
        previous = current
    end
    return area / 2, SVector(moment_x / 6, moment_y / 6)
end

function moments(c::StaticNgon)
    fallback = c.nr_verts == 0 ? (0u"m^2", SVector(0u"m^3", 0u"m^3")) :
        (zero(c.vertices[1].coords.x * c.vertices[1].coords.y),
         SVector(zero(c.vertices[1].coords.x^2 * c.vertices[1].coords.y),
             zero(c.vertices[1].coords.x * c.vertices[1].coords.y^2)))
    return _polygon_moments(c.vertices, c.nr_verts, fallback)
end

function _parabolic_correction(p::Parabola, v1::Point, v2::Point)
    dx1 = v1.coords.x - p.origin.coords.x
    dy1 = v1.coords.y - p.origin.coords.y
    dx2 = v2.coords.x - p.origin.coords.x
    dy2 = v2.coords.y - p.origin.coords.y
    τ1 = -p.𝛈[2] * dx1 + p.𝛈[1] * dy1
    τ2 = -p.𝛈[2] * dx2 + p.𝛈[1] * dy2
    Δτ = τ2 - τ1
    Δτ == zero(Δτ) && return zero(dx1 * dy1), SVector(
        zero(dx1 * dx1 * dy1), zero(dx1 * dy1 * dy1))

    η1 = p.𝛈[1] * dx1 + p.𝛈[2] * dy1 - p.shift
    η2 = p.𝛈[1] * dx2 + p.𝛈[2] * dy2 - p.shift
    a = (η2 - η1) / Δτ
    b = (η1 + η2) / 2 - a * (τ1 + τ2) / 2
    i0 = τ2 - τ1
    i1 = (τ2^2 - τ1^2) / 2
    i2 = (τ2^3 - τ1^3) / 3
    i3 = (τ2^4 - τ1^4) / 4
    i4 = (τ2^5 - τ1^5) / 5

    # Integral of (parabola - chord) in local (τ, η) coordinates.
    Δarea = -p.curvature * i2 / 2 - a * i1 - b * i0
    Δτmoment = -p.curvature * i3 / 2 - a * i2 - b * i1
    Δηmoment = p.curvature^2 * i4 / 8 -
        (a^2 * i2 + 2a * b * i1 + b^2 * i0) / 2
    Δnormal = p.shift * Δarea + Δηmoment
    Δx = p.origin.coords.x * Δarea + p.𝛈[1] * Δnormal - p.𝛈[2] * Δτmoment
    Δy = p.origin.coords.y * Δarea + p.𝛈[2] * Δnormal + p.𝛈[1] * Δτmoment
    return Δarea, SVector(Δx, Δy)
end

function moments(poly::StaticParabolicNgon)
    raw_area, raw_moment = _polygon_moments(poly.vertices, poly.nr_verts, poly.original_moments)
    if poly.nr_verts ≥ 2
        for index in 1:poly.nr_verts
            poly.parabolic_faces[index] || continue
            Δarea, Δmoment = _parabolic_correction(poly.curve, poly.vertices[index],
                poly.vertices[mod1(index + 1, poly.nr_verts)])
            raw_area += Δarea
            raw_moment += Δmoment
        end
    end
    return poly.complement ?
        (poly.original_moments[1] - raw_area, poly.original_moments[2] - raw_moment) :
        (raw_area, raw_moment)
end

smeasure(poly::StaticParabolicNgon) = moments(poly)[1]

@inline function _same_point(v1::Point, v2::Point; tol::Real=32eps())
    return abs(ustrip(v1.coords.x - v2.coords.x)) ≤ tol &&
           abs(ustrip(v1.coords.y - v2.coords.y)) ≤ tol
end

function _append_parabolic_vertex!(out::StaticParabolicNgon, vertex::Point,
    previous_is_parabolic::Bool)
    if out.nr_verts > 0 && _same_point(out.vertices[1], vertex)
        out.parabolic_faces[out.nr_verts] = previous_is_parabolic
        return out
    elseif out.nr_verts > 0 && _same_point(out.vertices[out.nr_verts], vertex)
        return out
    end
    out.nr_verts < capacity(out) || throw(ArgumentError(
        "StaticParabolicNgon capacity $(capacity(out)) is too small for this clipping operation"))
    if out.nr_verts > 0
        out.parabolic_faces[out.nr_verts] = previous_is_parabolic
    end
    out.nr_verts += 1
    out.vertices[out.nr_verts] = vertex
    out.parabolic_faces[out.nr_verts] = false
    return out
end

function _edge_roots(p::Parabola, v1::Point, v2::Point)
    dx = v2.coords.x - v1.coords.x
    dy = v2.coords.y - v1.coords.y
    rx = v1.coords.x - p.origin.coords.x
    ry = v1.coords.y - p.origin.coords.y
    τ0 = -p.𝛈[2] * rx + p.𝛈[1] * ry
    τd = -p.𝛈[2] * dx + p.𝛈[1] * dy
    A = p.curvature * τd^2 / 2
    B = p.𝛈[1] * dx + p.𝛈[2] * dy + p.curvature * τ0 * τd
    C = p.𝛈[1] * rx + p.𝛈[2] * ry - p.shift + p.curvature * τ0^2 / 2
    root1 = 0.0
    root2 = 0.0
    nroots = 0
    if A == zero(A)
        if B != zero(B)
            root1 = ustrip(-C / B)
            nroots = 1
        end
    else
        # Normalize before solving: the stable q-form avoids cancellation for
        # an interface almost parallel to an input edge.
        b = ustrip(B / A)
        c = ustrip(C / A)
        discriminant = b^2 - 4c
        if discriminant ≥ 0
            q = -(b + copysign(sqrt(discriminant), b)) / 2
            if q == 0
                root1 = -b / 2
                root2 = root1
            else
                root1 = q
                root2 = c / q
            end
            nroots = 2
        end
    end
    valid1 = nroots ≥ 1 && 32eps() < root1 < 1 - 32eps()
    valid2 = nroots == 2 && 32eps() < root2 < 1 - 32eps()
    if valid1 && valid2
        root1 > root2 && ((root1, root2) = (root2, root1))
        return abs(root1 - root2) ≤ 64eps() ? (root1, 0.0, 1) : (root1, root2, 2)
    elseif valid1
        return root1, 0.0, 1
    elseif valid2
        return root2, 0.0, 1
    end
    return 0.0, 0.0, 0
end

@inline _edge_point(v1::Point, v2::Point, s::Real) = v1 + s * (v2 - v1)

function _clip_negative_parabola!(out::StaticParabolicNgon, verts, p::Parabola)
    nverts = length(verts)
    _require_capacity(out, 2nverts + 2)
    out.nr_verts = 0
    nverts == 0 && return out
    all_inside = true
    start = nothing
    for index in eachindex(verts)
        value = distance(p, verts[index])
        all_inside &= value ≤ zero(value)
        value < zero(value) && isnothing(start) && (start = index)
    end
    if all_inside
        for vertex in verts
            _append_parabolic_vertex!(out, vertex, false)
        end
        out.nr_verts > 0 && (out.parabolic_faces[out.nr_verts] = false)
        return out
    end
    isnothing(start) && return out

    _append_parabolic_vertex!(out, verts[start], false)
    for offset in 0:nverts-1
        index = mod1(start + offset, nverts)
        next_index = mod1(index + 1, nverts)
        v1, v2 = verts[index], verts[next_index]
        inside = distance(p, v1) < zero(distance(p, v1))
        root1, root2, nroots = _edge_roots(p, v1, v2)
        for root_index in 1:nroots
            root = root_index == 1 ? root1 : root2
            next_root = root_index == nroots ? 1.0 : root2
            right_inside = distance(p, _edge_point(v1, v2, (root + next_root) / 2)) ≤
                zero(distance(p, v1))
            if inside != right_inside
                _append_parabolic_vertex!(out, _edge_point(v1, v2, root), !inside)
            end
            inside = right_inside
        end
        inside && _append_parabolic_vertex!(out, v2, false)
    end
    return out
end

"""
    intersect!(out, polygon, parabola)

Clip a polygon to the liquid side of a parabola.  The output retains which
faces are arcs, so `smeasure(out)` and `moments(out)` are exact up to floating
point round-off.  The implementation clips the complement for positive
curvature, where the retained side can otherwise be non-convex.
"""
function Base.intersect!(out::StaticParabolicNgon, c::Ngon, p::Parabola)
    out.original_moments = moments(c)
    out.curve = p
    out.complement = p.curvature > zero(p.curvature)
    pwork = out.complement ? complement(p) : p
    out.curve = pwork
    return _clip_negative_parabola!(out, c.vertices, pwork)
end

function Base.intersect!(out::StaticParabolicNgon, c::StaticNgon, p::Parabola)
    out.original_moments = moments(c)
    out.curve = p
    out.complement = p.curvature > zero(p.curvature)
    pwork = out.complement ? complement(p) : p
    out.curve = pwork
    return _clip_negative_parabola!(out, view(c.vertices, 1:c.nr_verts), pwork)
end

function Base.intersect(c::Ngon, p::Parabola)
    out = StaticParabolicNgon(c, p)
    return intersect!(out, c, p)
end

function smeasure(p::Parabola, c::Ngon;
    workspace::Union{Nothing,StaticParabolicNgon}=nothing)
    out = isnothing(workspace) ? StaticParabolicNgon(c, p) : workspace
    intersect!(out, c, p)
    return smeasure(out)
end

function smeasure(p::Parabola, c::StaticNgon;
    workspace::Union{Nothing,StaticParabolicNgon}=nothing)
    out = isnothing(workspace) ? StaticParabolicNgon(c, p) : workspace
    intersect!(out, c, p)
    return smeasure(out)
end

function moments(p::Union{PlanarHS{2},Parabola}, c::Ngon;
    workspace=nothing)
    if p isa PlanarHS
        out = isnothing(workspace) ? StaticNgon(c) : workspace
        intersect!(out, c, p)
        return moments(out)
    end
    out = isnothing(workspace) ? StaticParabolicNgon(c, p) : workspace
    intersect!(out, c, p)
    return moments(out)
end

function _parabolic_shift(c::Ngon, 𝛈::SVector{2}, curvature::Quantity,
    origin::Point, target_area::Quantity; workspace::Union{Nothing,StaticParabolicNgon}=nothing,
    iterations::Integer=56)
    total_area = smeasure(c)
    zero(total_area) ≤ target_area ≤ total_area || throw(DomainError(target_area,
        "area must lie between zero and the polygon's signed area"))
    τmax = zero(c.vertices[1].coords.x - origin.coords.x)
    ηmin = 𝛈[1] * (c.vertices[1].coords.x - origin.coords.x) +
           𝛈[2] * (c.vertices[1].coords.y - origin.coords.y)
    ηmax = ηmin
    for vertex in c.vertices
        dx, dy = vertex.coords.x - origin.coords.x, vertex.coords.y - origin.coords.y
        τ = -𝛈[2] * dx + 𝛈[1] * dy
        η = 𝛈[1] * dx + 𝛈[2] * dy
        abs(τ) > τmax && (τmax = abs(τ))
        η < ηmin && (ηmin = η)
        η > ηmax && (ηmax = η)
    end
    padding = abs(curvature) * τmax^2 / 2
    lower, upper = ηmin - padding, ηmax + padding
    target_area == zero(target_area) && return lower
    target_area == total_area && return upper
    p = Parabola(𝛈, lower, curvature, origin)
    out = isnothing(workspace) ? StaticParabolicNgon(c, p) : workspace
    for _ in 1:iterations
        middle = (lower + upper) / 2
        p = Parabola(𝛈, middle, curvature, origin)
        area = smeasure(p, c; workspace=out)
        if area < target_area
            lower = middle
        else
            upper = middle
        end
    end
    return (lower + upper) / 2
end

function _parabola_with_area(θ::Real, curvature::Quantity, area::Quantity, c::Ngon,
    origin::Point; workspace::Union{Nothing,StaticParabolicNgon}=nothing)
    𝛈 = angle_to_normal(θ)
    shift = _parabolic_shift(c, 𝛈, curvature, origin, area; workspace=workspace)
    return Parabola(𝛈, shift, curvature, origin)
end
