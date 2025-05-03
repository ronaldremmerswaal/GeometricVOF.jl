function reconstruct(p0::PlanarHS{2}, α_central::T, c_central::Ngon, αs::AbstractArray{T}, cs::SubDomain,
    cmeasures::AbstractArray{Q}=measure.(cs); workspace::StaticNgon=StaticNgon(c_central)) where {T <: Real, Q <: Quantity}
    ref_vol = smeasure(c_central) * α_central

    wrapped_costfun(θ::Real) = lvira_costfun(PlanarHS(θ, ref_vol, c_central; workspace=workspace), cs, αs, cmeasures, c_central, workspace=workspace)

    θ0 = GeometricVOF.normal_to_angle(p0.𝛈)
    θ = brent_min(wrapped_costfun, θ0; xtol=1E-8, maxiters=25, step_max=.5)

    return PlanarHS(θ, ref_vol, c_central; workspace=workspace)
end

function lvira_costfun(p::PlanarHS{2}, cs::SubDomain, αs::AbstractArray{T}, cmeasures::AbstractArray{Q}, c_central::Ngon; workspace::StaticNgon=StaticNgon(c_central)) where {T <: Real, Q <: Quantity}
    err = 0    # Value
    derr = 0   # Derivative w.r.t. the angle

    intersect!(workspace, c_central, p)

    c_iface = Segment(workspace.vertices[workspace.interface_index], workspace.vertices[mod1(workspace.interface_index + 1, workspace.nr_verts)]) # NOTE: this assumes the polygon is convex
    # c_iface_area = measure(c_iface)
    c_iface_centroid = centroid(c_iface)
    tangent = [-p.𝛈[2], p.𝛈[1]]

    dshift = tangent ⋅ to(c_iface_centroid) # The shift derivative that ensures that the central volume is invariant
    for (c, α, cmeas) ∈ zip(cs, αs, cmeasures)
        intersect!(workspace, c, p)
        if workspace.nr_verts < 3
            continue
        end
        err_local = smeasure(workspace) / cmeas - α

        # If workspace.interface_index == 0 then there is no interface inside this cell
        if workspace.interface_index > 0
            iface = Segment(workspace.vertices[workspace.interface_index], workspace.vertices[mod1(workspace.interface_index + 1, workspace.nr_verts)]) # NOTE: this assumes the polygon is convex
            iface_area = measure(iface)
            iface_centroid = centroid(iface)

            derr_local = iface_area * (dshift - tangent ⋅ to(iface_centroid))
            derr += + 2 * err_local * derr_local / cmeas
        end

        err += err_local^2
    end

    return err, derr
end
