function reconstruct(cfun::CostFunction, c::Ngon, α::Quantity, p0::PlanarHS{2})
    function wrapped_cfun(θ::Vector)
        𝛈 = GeometricVOF.angle_to_normal(θ[1])
        s = shift(c, 𝛈, α)
        cfun(PlanarHS{2}(𝛈, s))
    end

    θ0 = GeometricVOF.normal_to_angle(p0.𝛈)
    res = optimize(wrapped_cfun, [θ0], g_tol=1E-30, x_tol=1E-6)
    # TODO implement linearization of cost function

    # TODO make a conversion function: (θ, α) -> PlanarHS
    θ = res.minimizer[1]
    𝛈 = GeometricVOF.angle_to_normal(θ)
    s = shift(c, 𝛈, α)

    return PlanarHS{2}(𝛈, s)
end

abstract type CostFunction end

struct LVIRA <: CostFunction
    cs::Vector{Ngon}
    αs::Vector{Quantity}
    scalings::Vector{Quantity}
end
LVIRA(cs::Vector{N}, αs::Vector{T}) where {N <: Ngon, T <: Quantity} =
    LVIRA(cs, αs, measure.(cs))

function (f::LVIRA)(p::PlanarHS{2})
    err0 = 0
    for (c, α, scaling) ∈ zip(f.cs, f.αs, f.scalings)
        err0 += ((measure(p, c) - α) / scaling)^2
    end
    return err0
end
