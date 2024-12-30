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
        M += measure(Φ, c_not_p)
    end

    if !isnothing(c_p)
        M += measure(not_Φ, c_p)
    end

    return M
end
