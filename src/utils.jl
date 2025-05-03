"""
    symmetric_difference(Φ, p, c)

Computes the area of the symmetric difference between the halfspace defined by Φ and p
within c. I.e. the following symmetric difference is computed
    {𝐱 ∈ c | Φ(𝐱) ≤ 0} Δ {𝐱 ∈ c | p.𝛈 ⋅ 𝐱 - p.shift ≤ 0}
"""
function symmetric_difference(Φ::Function, p::PlanarHS{2}, c::Ngon)
    c_p = c ∩ p
    c_not_p = c ∩ complement(p)

    not_Φ(x, y) = -Φ(x, y)  # TODO doesn't work for Bool

    M = 0u"m^2"
    if !isnothing(c_not_p)
        M += smeasure(Φ, c_not_p)
    end

    if !isnothing(c_p)
        M += smeasure(not_Φ, c_p)
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
    mom0 = 0u"m^2"
    mom1 = [0, 0]u"m^3"
    verts = vertices(c)

    for vdx ∈ eachindex(verts)
        v1 = verts[vdx]
        v2 = verts[vdx==length(verts) ? 1 : vdx + 1]
        mom0_tmp = v1.coords.x * v2.coords.y - v1.coords.y * v2.coords.x

        mom0 += mom0_tmp
        mom1[1] += mom0_tmp * (v1.coords.x + v2.coords.x)
        mom1[2] += mom0_tmp * (v1.coords.y + v2.coords.y)
    end

    return mom0 / 2, mom1 / 6
end

angle_to_normal(θ::Number) = SVector{2}(cos(θ), sin(θ))
normal_to_angle(𝛈::AbstractVector) = atan(𝛈[2], 𝛈[1])
