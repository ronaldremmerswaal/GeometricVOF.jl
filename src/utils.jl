"""
    symmetric_difference(Φ, p, c)

Computes the area of the symmetric difference between the halfspace defined by Φ and p
within c. I.e. the following symmetric difference is computed
    {𝐱 ∈ c | Φ(𝐱) ≤ 0} Δ {𝐱 ∈ c | p.𝛈 ⋅ 𝐱 - p.shift ≤ 0}
"""
function symmetric_difference(Φ::F, p::PlanarHS{2}, c::Ngon) where {F}
    c_p = c ∩ p
    c_not_p = c ∩ complement(p)

    not_Φ(x, y) = -Φ(x, y)  # TODO doesn't work for Bool

    M = 0u"m^2"
    if !isnothing(c_not_p)
        M += abs(smeasure(Φ, c_not_p))
    end

    if !isnothing(c_p)
        M += abs(smeasure(not_Φ, c_p))
    end

    return M
end

"""
    symmetric_difference(p1, p2, c)

Computes the area of the symmetric difference between the halfspace defined by p1 and p2
within c. I.e. the following symmetric difference is computed
    {𝐱 ∈ c | p1.𝛈 ⋅ 𝐱 - p1.shift ≤ 0} Δ {𝐱 ∈ c | p2.𝛈 ⋅ 𝐱 - p2.shift ≤ 0}
"""
function symmetric_difference(p1::PlanarHS{2}, p2::PlanarHS{2}, c::Ngon)
    c_p1 = c ∩ p1
    c_not_p1 = c ∩ complement(p1)

    M = 0u"m^2"
    if !isnothing(c_not_p1)
        M += abs(smeasure(p2, c_not_p1))
    end

    if !isnothing(c_p1)
        M += abs(smeasure(complement(p2), c_p1))
    end

    return M
end

"""
    moments(c)

Zeroth and first order moments of Ngon `c`.

# Examples
```julia-repl
julia> poly = Triangle((0, 0), (1, 0), (0, 1))
julia> moments(poly)
(0.5 m^2, [0.16666666666666666 m^3, 0.16666666666666666 m^3])
```
"""
function moments(c::Ngon)
    verts = vertices(c)
    fallback = (zero(verts[1].coords.x * verts[1].coords.y),
        SVector(zero(verts[1].coords.x^2 * verts[1].coords.y),
            zero(verts[1].coords.x * verts[1].coords.y^2)))
    return _polygon_moments(verts, length(verts), fallback)
end

angle_to_normal(θ::Number) = SVector{2}(cos(θ), sin(θ))
normal_to_angle(𝛈::AbstractVector) = atan(𝛈[2], 𝛈[1])
