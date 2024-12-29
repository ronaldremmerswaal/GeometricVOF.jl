function reconstruct(cfun::Function, c::Ngon, α::Quantity, p0::PlanarHS{2})
    function wrapped_cfun(θ)
        𝛈 = GeometricVOF.angle_to_normal(θ)
        s = shift(c, 𝛈, α)
        cfun(PlanarHS{2}(𝛈, s))
    end


end

abstract type CostFunction end

struct LVIRA <: CostFunction
    cs::Vector{Ngon}
    αs::Vector{Quantity}
    scalings::Vector{Quantity}
end
LVIRA(cs::Vector{Ngon}, αs::Vector) = LVIRA(cs, αs, measure.(cs))

function (f::LVIRA)(p::PlanarHS{2})
    err0 = 0
    for (c, α, scaling) ∈ zip(f.cs, f.αs, f.scalings)
        err0 += ((measure(c, p) - α) / scaling)^2
    end
    return err0
end
