"""
    symmetric_difference(Φ, p, c)

Computes the area of the symmetric difference between the halfspace defined by Φ and p
within c. I.e. the following symmetric difference is computed
    {𝐱 ∈ c | Φ(𝐱) ≤ 0} Δ {𝐱 ∈ c | p.𝛈 ⋅ 𝐱 - p.shift ≤ 0}
"""
function symmetric_difference(Φ::Function, p::PlanarHS{2}, c::Ngon)
    c_p = c ∩ p
    c_not_p = c ∩ complement(p)

    not_Φ(x, y) = -Φ(x, y)

    measure(Φ, c_not_p) + measure(not_Φ, c_p)
end
