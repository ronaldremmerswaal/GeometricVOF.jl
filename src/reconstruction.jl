function reconstruct(p0::PlanarHS{2}, α_central::T, c_central::Ngon, αs::AbstractArray{T}, cs::SubDomain,
    cmeasures::AbstractArray{Q}=smeasure.(cs); verbose::Bool=false, xatol::Real=√(eps(T)), workspace::StaticNgon=StaticNgon(c_central), shift_workspace::MVector=MVector{32, Float64}(undef)) where {T <: Real, Q <: Quantity}
    ref_vol = smeasure(c_central) * α_central

    wrapped_costfun(θ::Real) = lvira_costfun(PlanarHS(θ, ref_vol, c_central; workspace=workspace, shift_workspace=shift_workspace), cs, αs, cmeasures, c_central, workspace=workspace)

    θ0 = GeometricVOF.normal_to_angle(p0.𝛈)
    θ = brent_min(wrapped_costfun, θ0; xatol=xatol, maxiters=25, step_max=.5, verbose=verbose)

    return PlanarHS(θ, ref_vol, c_central; workspace=workspace, shift_workspace=shift_workspace)
end

function lvira_costfun(p::PlanarHS{2}, cs::SubDomain, αs::AbstractArray{T}, cmeasures::AbstractArray{Q}, c_central::Ngon; workspace::StaticNgon=StaticNgon(c_central)) where {T <: Real, Q <: Quantity}
    err  = 0.0   # initialise as Float64 to avoid Int→Float type change mid-loop
    derr = 0.0

    intersect!(workspace, c_central, p)

    # Centroid of the central interface edge — computed inline to avoid Segment/centroid allocs
    ci1 = workspace.vertices[workspace.interface_index]
    ci2 = workspace.vertices[mod1(workspace.interface_index + 1, workspace.nr_verts)]
    tangent = SVector(-p.𝛈[2], p.𝛈[1])

    dshift = tangent[1] * (ci1.coords.x + ci2.coords.x) / 2 +
             tangent[2] * (ci1.coords.y + ci2.coords.y) / 2  # ≡ tangent ⋅ centroid(c_iface)
    for (c, α, cmeas) ∈ zip(cs, αs, cmeasures)
        intersect!(workspace, c, p)
        if workspace.nr_verts < 3
            continue
        end
        err_local = smeasure(workspace) / cmeas - α

        ω = 1 / (α * (1 - α) + 1E-2)

        # If workspace.interface_index == 0 then there is no interface inside this cell
        if workspace.interface_index > 0
            # Inline segment measure and centroid — avoids Segment + centroid + to() allocations
            iv1 = workspace.vertices[workspace.interface_index]
            iv2 = workspace.vertices[mod1(workspace.interface_index + 1, workspace.nr_verts)]
            idx = iv2.coords.x - iv1.coords.x
            idy = iv2.coords.y - iv1.coords.y
            iface_area = sqrt(idx^2 + idy^2)   # ≡ Meshes.measure(Segment(iv1,iv2))
            iface_cx   = (iv1.coords.x + iv2.coords.x) / 2
            iface_cy   = (iv1.coords.y + iv2.coords.y) / 2

            derr_local = iface_area * (dshift - (tangent[1] * iface_cx + tangent[2] * iface_cy))
            derr += 2ω * err_local * derr_local / cmeas
        end

        err += ω * err_local^2
    end

    return err, derr
end
