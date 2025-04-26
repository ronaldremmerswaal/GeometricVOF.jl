abstract type CostFunction end

function reconstruct(cfun::CostFunction, c::Ngon, αvol::Quantity, p0::PlanarHS{2})
    wrapped_cfun(θ::Real) = cfun(PlanarHS(θ, αvol, c))

    θ0 = GeometricVOF.normal_to_angle(p0.𝛈)
    θ = brent_min(wrapped_cfun, θ0; xtol=1E-8, maxiters=25, step_max=.5)

    return PlanarHS(θ, αvol, c)
end

struct LVIRA{T} <: CostFunction
    cs::AbstractArray{Ngon}
    αs::AbstractArray{T}
    cmeasures::AbstractArray{Quantity}

    c_c::Ngon # The `central' cell`
end
LVIRA(cs::AbstractArray{N}, αs::AbstractArray{T}, c_c::Ngon) where {N <: Ngon, T <: Real} =
    LVIRA{T}(cs, αs, measure.(cs), c_c)

function (f::LVIRA)(p::PlanarHS{2})
    err0 = 0    # Value
    derr0 = 0   # Derivative w.r.t. the angle

    c_cp = f.c_c ∩ p
    c_iface = Segment(c_cp.vertices[1], c_cp.vertices[2]) # NOTE: this assumes the polygon is convex
    # c_iface_area = measure(c_iface)
    c_iface_centroid = centroid(c_iface)
    c_iface_tangent = tangent(c_iface)

    dshift = c_iface_tangent ⋅ to(c_iface_centroid) # The shift derivative that ensures that the central volume is invariant
    for (c, α, cmeas) ∈ zip(f.cs, f.αs, f.cmeasures)
        cp = c ∩ p
        if isnothing(cp)
            continue
        end
        err_local = measure(cp) / cmeas - α

        iface = Segment(cp.vertices[1], cp.vertices[2]) # NOTE: this assumes the polygon is convex
        iface_area = measure(iface)
        iface_centroid = centroid(iface)

        derr0_local = iface_area * (dshift - c_iface_tangent ⋅ to(iface_centroid))

        err0 += err_local^2
        derr0 = derr0 + 2 * err_local * derr0_local / cmeas
    end
    return err0, derr0
end
