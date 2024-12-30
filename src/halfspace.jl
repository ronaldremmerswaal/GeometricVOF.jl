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
    first_inside = 0
    for (vdx, vert) ∈ enumerate(c.vertices)
        dist[vdx] = distance(p, vert)
        if dist[vdx] ≤ 0u"m"
            nr_inside += 1
            if first_inside == 0
                first_inside = vdx
            end
        end
    end

    # Quick return in trivial cases
    if nr_inside == 0
        return nothing  # TODO type stability
    elseif nr_inside == length(c.vertices)
        return c        # TODO no copy available
    end

    new_verts = Vector{eltype(c.vertices)}()

    # Construct new polygon by looping over the edges of the old polygon
    vdx = first_inside
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
    shift0 = 𝛈 ⋅ to(centroid(c))
    α_err0 = α_err(shift0)

    if α_err0 == 0u"m" return shift0 end

    shift_min, shift_max = shift_extrema(c, 𝛈)

    bracket = α_err0 > 0u"m^2" ? (shift_min, shift0) : (shift0, shift_max)

    return find_zero(α_err, bracket, Roots.Brent())
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

angle_to_normal(θ::Number) = SVector{2}(cos(θ), sin(θ))
normal_to_angle(𝛈::AbstractVector) = atan(𝛈[2], 𝛈[1])
