function reconstruct(p0::PlanarHS{2}, α_central::T, c_central::Ngon, cs::SubDomain, αs::AbstractArray{T}, cmeasures::AbstractArray{Q}=measure.(cs)) where {T <: Real, Q <: Quantity}
    ref_vol = measure(c_central) * α_central
    wrapped_costfun(θ::Real) = lvira_costfun(PlanarHS(θ, ref_vol, c_central), cs, αs, cmeasures, c_central)

    θ0 = GeometricVOF.normal_to_angle(p0.𝛈)
    θ = brent_min(wrapped_costfun, θ0; xtol=1E-8, maxiters=25, step_max=.5)

    return PlanarHS(θ, ref_vol, c_central)
end

function lvira_costfun(p::PlanarHS{2}, cs::SubDomain, αs::AbstractArray{T}, cmeasures::AbstractArray{Q}, c_central::Ngon) where {T <: Real, Q <: Quantity}
    err = 0    # Value
    derr = 0   # Derivative w.r.t. the angle

    cp_memory = MVector{30, eltype(c_central.vertices)}(undef)
    cp_memory, cp_length, cp_interface = intersect!(cp_memory, c_central, p)

    c_iface = Segment(cp_memory[cp_interface], cp_memory[mod1(cp_interface + 1, cp_length)]) # NOTE: this assumes the polygon is convex
    # c_iface_area = measure(c_iface)
    c_iface_centroid = centroid(c_iface)
    tangent = [-p.𝛈[2], p.𝛈[1]]

    dshift = tangent ⋅ to(c_iface_centroid) # The shift derivative that ensures that the central volume is invariant
    for (c, α, cmeas) ∈ zip(cs, αs, cmeasures)
        cp_memory, cp_length, cp_interface = intersect!(cp_memory, c, p)
        if cp_length < 3
            continue
        end
        err_local = measure(cp_memory, cp_length) / cmeas - α

        # If cp_interface == 0 then there is no interface inside this cell
        if cp_interface > 0
            iface = Segment(cp_memory[cp_interface], cp_memory[mod1(cp_interface + 1, cp_length)]) # NOTE: this assumes the polygon is convex
            iface_area = measure(iface)
            iface_centroid = centroid(iface)

            derr_local = iface_area * (dshift - tangent ⋅ to(iface_centroid))
            derr += + 2 * err_local * derr_local / cmeas
        end

        err += err_local^2
    end

    return err, derr
end
