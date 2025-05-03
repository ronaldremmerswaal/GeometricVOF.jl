abstract type CostFunction end

function reconstruct(cfun::CostFunction, αvol::Quantity, p0::PlanarHS{2})
    wrapped_cfun(θ::Real) = cfun(PlanarHS(θ, αvol, cfun.c_c))

    θ0 = GeometricVOF.normal_to_angle(p0.𝛈)
    θ = brent_min(wrapped_cfun, θ0; xtol=1E-8, maxiters=25, step_max=.5)

    return PlanarHS(θ, αvol, cfun.c_c)
end

struct LVIRA{T} <: CostFunction
    cs::SubDomain
    αs::AbstractArray{T}
    cmeasures::AbstractArray{Quantity}

    c_c::Ngon # The `central' cell`
end
LVIRA(cs::SubDomain, αs::AbstractArray{T}, c_c::Ngon) where {T <: Real} =
    LVIRA{T}(cs, αs, measure.(cs), c_c)

function (f::LVIRA)(p::PlanarHS{2})
    err = 0    # Value
    derr = 0   # Derivative w.r.t. the angle

    cp_memory = MVector{30, eltype(f.c_c.vertices)}(undef)
    (cp_memory, cp_length) = intersect!(cp_memory, f.c_c, p)

    # c_iface = Segment(c_cp.vertices[1], c_cp.vertices[2]) # NOTE: this assumes the polygon is convex
    # # c_iface_area = measure(c_iface)
    # c_iface_centroid = centroid(c_iface)
    # tangent = [-p.𝛈[2], p.𝛈[1]]


    # dshift = tangent ⋅ to(c_iface_centroid) # The shift derivative that ensures that the central volume is invariant
    for (c, α, cmeas) ∈ zip(f.cs, f.αs, f.cmeasures)
        (cp_memory, cp_length) = intersect!(cp_memory, c, p)
        if cp_length == 0
            continue
        end
        err_local = measure(cp_memory, cp_length) / cmeas - α

        # TODO first edge no longer interface
        # iface = Segment(cp_memory[1], cp_memory[2]) # NOTE: this assumes the polygon is convex
        # iface_area = measure(iface)
        # iface_centroid = centroid(iface)

        # derr_local = iface_area * (dshift - tangent ⋅ to(iface_centroid))

        err += err_local^2
        # derr += + 2 * err_local * derr_local / cmeas
    end
    return err, derr
end
