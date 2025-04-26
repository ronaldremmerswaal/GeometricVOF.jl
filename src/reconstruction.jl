abstract type CostFunction end

function reconstruct(cfun::CostFunction, c::Ngon, αvol::Quantity, p0::PlanarHS{2})
    wrapped_cfun(θ::Vector) = cfun(PlanarHS(θ[1], αvol, c))

    θ0 = GeometricVOF.normal_to_angle(p0.𝛈)
    res = optimize(wrapped_cfun, [θ0], g_tol=1E-30, x_abstol=1E-6)
    # TODO implement linearization of cost function

    # TODO make a conversion function: (θ, α) -> PlanarHS
    θ = res.minimizer[1]

    return PlanarHS(θ, αvol, c)
end

struct LVIRA{T} <: CostFunction
    cs::AbstractArray{Ngon}
    αs::AbstractArray{T}
    cmeasures::AbstractArray{Quantity}
end
LVIRA(cs::AbstractArray{N}, αs::AbstractArray{T}) where {N <: Ngon, T <: Real} =
    LVIRA{T}(cs, αs, measure.(cs))

function (f::LVIRA)(p::PlanarHS{2})
    err0 = 0
    for (c, α, cmeas) ∈ zip(f.cs, f.αs, f.cmeasures)
        err0 += (measure(p, c) / cmeas - α)^2
    end
    return err0
end
