abstract type HalfSpace{D} end

"""
    Plane(𝛈, shift)

Represents the half-space defined by

    𝛈 ⋅ 𝐱 ≤ shift
"""
struct Plane{D} <: HalfSpace{D}
    𝛈::SVector{D}
    shift::Number
end
Plane(𝛈::Vector, shift::Number) = Plane{length(𝛈)}(SVector{length(𝛈)}(𝛈), shift)

distance(p::Plane, 𝐱::Point) = p.𝛈 ⋅ 𝐱 - p.shift


"""
    intersect!(c, plane)

Intersects `c` with halfspace defined by `plane`.

# Examples
```julia-repl
julia> c = Triangle((0., 0.), (1., 0.), (0., 1.))
julia> p = Plane([-1., 0.], -0.5)
julia> intersect(c, p)
Triangle((0.5, 0.5), (0.5, 0.0), (1.0, 0.0)) # TODO check
```
"""
function Base.intersect(c::Ngon, p::Plane{2})

    nr_old_verts = length(vertices(c))

    # Compute the distance of each vertex to the plane
    dist = zeros(nr_old_verts)u"m"
    nr_inside = 0
    first_inside = 0
    for (vdx, vert) ∈ enumerate(c.vertices)
        dist[vdx] = distance(p, vert)
        if dist[vdx] ≤ 0
            nr_inside += 1
            if first_inside == 0
                first_inside = vdx
            end
        end
    end

    # Quick return in trivial cases
    if nr_inside == 0
        return nothing # TODO type stability
    elseif nr_inside == length(c.vertices)
        return copy(c)
    end

    new_verts = Vector{eltype(c.vertices)}()

    # Construct new polygon by looping over the edges of the old polygon
    vdx = first_inside
    for _ ∈ 1:nr_old_verts
        ndx = vdx == nr_old_verts ? 1 : vdx + 1
        edge_is_bisected = (dist[vdx] ≤ 0u"m") != (dist[ndx] ≤ 0u"m")

        if edge_is_bisected
            coeff = abs(dist[vdx] / (dist[ndx] - dist[vdx]))
            if (coeff < 1 || dist[ndx] > 0) && (coeff > 0 || dist[vdx] > 0)
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
        return Ngon(new_verts)
    end
end
