"""
    lvira(initial, fraction, cell, neighbor_fractions, neighbor_cells;
          cell_areas=smeasure.(neighbor_cells), workspace, shift_workspace)

Reconstruct a planar interface in `cell` with LVIRA. `fraction` must be
strictly between zero and one. Reuse the optional workspaces in a cell loop.
"""
function lvira(p0::PlanarHS{2}, α_central::T, c_central::Ngon,
    αs::AbstractArray{T}, cs::Union{SubDomain, AbstractArray{<:Ngon}},
    cmeasures::AbstractArray{Q}=smeasure.(cs);
    verbose::Bool=false,
    xatol::Real=√eps(T),
    workspace::StaticNgon=StaticNgon(c_central),
    shift_workspace::AbstractVector{<:Real}=Vector{Float64}(undef, length(c_central.vertices)),
) where {T<:Real, Q<:Quantity}
    _validate_reconstruction_inputs(α_central, αs, cs, cmeasures)
    ref_vol = smeasure(c_central) * α_central

    wrapped_costfun(θ::Real) = lvira_costfun(PlanarHS(θ, ref_vol, c_central; workspace=workspace, shift_workspace=shift_workspace), cs, αs, cmeasures, c_central, workspace=workspace)

    θ0 = GeometricVOF.normal_to_angle(p0.𝛈)
    θ = brent_min(wrapped_costfun, θ0; xatol=xatol, maxiters=25, step_max=.5, verbose=verbose)

    return PlanarHS(θ, ref_vol, c_central; workspace=workspace, shift_workspace=shift_workspace)
end

function _validate_reconstruction_inputs(α_central, αs, cs, cmeasures)
    zero(α_central) < α_central < one(α_central) || throw(DomainError(α_central,
        "LVIRA reconstruction requires a volume fraction strictly between zero and one"))
    length(αs) == length(cs) == length(cmeasures) || throw(DimensionMismatch(
        "neighbor fractions, cells, and cell areas must have the same length"))
    return nothing
end

function lvira_costfun(p::PlanarHS{2}, cs::Union{SubDomain, AbstractArray{<:Ngon}}, αs::AbstractArray{T},
    cmeasures::AbstractArray{Q}, c_central::Ngon;
    workspace::StaticNgon=StaticNgon(c_central),
) where {T<:Real, Q<:Quantity}
    _lvira_costfun(p, cs, αs, cmeasures, c_central; workspace=workspace)
end

function lvira_costfun(p::PlanarHS{2}, cs::AbstractVector{<:StaticNgon},
    αs::AbstractVector{T}, cmeasures::AbstractVector{Q}, c_central::StaticNgon;
    workspace::StaticNgon=StaticNgon(eltype(c_central.vertices), capacity(c_central) + 2),
) where {T<:Real, Q<:Quantity}
    _lvira_costfun(p, cs, αs, cmeasures, c_central; workspace=workspace)
end

function _lvira_costfun(p::PlanarHS{2}, cs, αs, cmeasures, c_central; workspace)
    err  = 0.0
    derr = 0.0

    intersect!(workspace, c_central, p)

    # Centroid of the central interface edge — computed inline to avoid Segment/centroid allocs
    ci1 = workspace.vertices[workspace.interface_index]
    ci2 = workspace.vertices[mod1(workspace.interface_index + 1, workspace.nr_verts)]
    tangent = SVector(-p.𝛈[2], p.𝛈[1])

    dshift = tangent[1] * (ci1.coords.x + ci2.coords.x) / 2 +
             tangent[2] * (ci1.coords.y + ci2.coords.y) / 2  # tangent ⋅ centroid(c_iface)
    for (c, α, cmeas) ∈ zip(cs, αs, cmeasures)
        intersect!(workspace, c, p)
        if workspace.nr_verts < 3
            continue
        end
        err_local = smeasure(workspace) / cmeas - α

        ω = 1 / (α * (1 - α) + 1E-2)

        # If workspace.interface_index == 0 then there is no interface inside this cell
        if workspace.interface_index > 0
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

@inline function _validate_fraction(α)
    zero(α) < α < one(α) || throw(DomainError(α,
        "interface reconstruction requires a volume fraction strictly between zero and one"))
    return nothing
end

function _parabolic_monomials(poly::StaticParabolicNgon)
    p = poly.curve
    τ0 = zero(p.shift)
    i0 = zero(τ0)
    i1 = zero(τ0^2)
    i2 = zero(τ0^3)
    i3 = zero(τ0^4)
    i4 = zero(τ0^5)
    i5 = zero(τ0^6)
    for index in 1:poly.nr_verts
        poly.parabolic_faces[index] || continue
        v1 = poly.vertices[index]
        v2 = poly.vertices[mod1(index + 1, poly.nr_verts)]
        τ1 = -p.𝛈[2] * (v1.coords.x - p.origin.coords.x) + p.𝛈[1] * (v1.coords.y - p.origin.coords.y)
        τ2 = -p.𝛈[2] * (v2.coords.x - p.origin.coords.x) + p.𝛈[1] * (v2.coords.y - p.origin.coords.y)
        i0 += τ2 - τ1
        i1 += (τ2^2 - τ1^2) / 2
        i2 += (τ2^3 - τ1^3) / 3
        i3 += (τ2^4 - τ1^4) / 4
        i4 += (τ2^5 - τ1^5) / 5
        i5 += (τ2^6 - τ1^6) / 6
    end
    return i0, i1, i2, i3, i4, i5
end

function _parabolic_shift_derivatives(poly::StaticParabolicNgon)
    i0, i1, i2, i3, _, _ = _parabolic_monomials(poly)
    p = poly.curve
    i0 == zero(i0) && return zero(p.shift), zero(p.shift^2)
    raw_angle = (i1 - i1 * p.shift * p.curvature + i3 * p.curvature^2 / 2) / i0
    curvature = i2 / (2i0)
    return poly.complement ? -raw_angle : raw_angle, curvature
end

function _parabolic_volume_derivatives(poly::StaticParabolicNgon, dshift_angle, dshift_curvature)
    i0, i1, i2, i3, _, _ = _parabolic_monomials(poly)
    p = poly.curve
    dshift = poly.complement ? -dshift_angle : dshift_angle
    raw_angle = i0 * dshift - i1 + i1 * p.shift * p.curvature - i3 * p.curvature^2 / 2
    curvature = i0 * dshift_curvature - i2 / 2
    return poly.complement ? -raw_angle : raw_angle, curvature
end

function _parabolic_first_moment_angle_derivative(poly::StaticParabolicNgon, dshift_angle)
    i0, i1, i2, i3, i4, i5 = _parabolic_monomials(poly)
    p = poly.curve
    dshift = poly.complement ? -dshift_angle : dshift_angle
    dη = p.shift * dshift * i0 - p.shift * i1 + i3 * p.curvature / 2 -
        i2 * dshift * p.curvature / 2 + p.curvature *
        (i1 * p.shift^2 + i5 * p.curvature^2 / 4 - p.shift * p.curvature * i3)
    dτ = dshift * i1 - i2 + p.curvature * p.shift * i2 - i4 * p.curvature^2 / 2
    derivative = SVector(p.𝛈[1] * dη - p.𝛈[2] * dτ,
        p.𝛈[2] * dη + p.𝛈[1] * dτ)
    return poly.complement ? -derivative : derivative
end

function _planar_first_moment_angle_derivative(poly::StaticNgon, p::PlanarHS)
    poly.interface_index == 0 && return SVector(zero(p.shift^3), zero(p.shift^3))
    v1 = poly.vertices[poly.interface_index]
    v2 = poly.vertices[mod1(poly.interface_index + 1, poly.nr_verts)]
    τ1 = -p.𝛈[2] * v1.coords.x + p.𝛈[1] * v1.coords.y
    τ2 = -p.𝛈[2] * v2.coords.x + p.𝛈[1] * v2.coords.y
    i0 = τ2 - τ1
    i0 == zero(i0) && return SVector(zero(p.shift^3), zero(p.shift^3))
    i1 = (τ2^2 - τ1^2) / 2
    i2 = (τ2^3 - τ1^3) / 3
    dshift = i1 / i0
    dη = p.shift * dshift * i0 - p.shift * i1
    dτ = dshift * i1 - i2
    return SVector(p.𝛈[1] * dη - p.𝛈[2] * dτ,
        p.𝛈[2] * dη + p.𝛈[1] * dτ)
end

"""
    mof(initial, fraction, first_moment, cell; workspace, shift_workspace)

Moment-of-fluid (MOF) reconstruction of a planar interface.  `first_moment`
is the liquid first moment (not its centroid) in the global coordinate system.
The returned `PlanarHS` has the requested liquid volume and minimizes the
Euclidean first-moment mismatch.
"""
function mof(p0::PlanarHS{2}, α::Real, first_moment, c::Ngon;
    workspace::StaticNgon=StaticNgon(c),
    shift_workspace::AbstractVector{<:Real}=Vector{Float64}(undef, length(c.vertices)),
    xatol::Real=√eps(Float64),
    maxiters::Integer=25,
    step_max::Real=.5,
    verbose::Bool=false,
)
    _validate_fraction(α)
    cell_area = abs(smeasure(c))
    target_area = α * cell_area
    moment_scale = cell_area * sqrt(cell_area)
    cost_and_derivative(θ) = begin
        p = PlanarHS(θ, target_area, c; workspace=workspace, shift_workspace=shift_workspace)
        intersect!(workspace, c, p)
        _, reconstructed_moment = moments(workspace)
        difference = (reconstructed_moment - first_moment) / moment_scale
        derivative = _planar_first_moment_angle_derivative(workspace, p) / moment_scale
        return sum(abs2, difference), 2 * dot(difference, derivative)
    end
    θ = brent_min(cost_and_derivative, normal_to_angle(p0.𝛈);
        xatol=xatol, maxiters=maxiters, step_max=step_max, verbose=verbose)
    return PlanarHS(θ, target_area, c; workspace=workspace, shift_workspace=shift_workspace)
end

function _parabolic_workspace(c::Ngon, p::Parabola, cs=())
    nvertices = length(c.vertices)
    for neighbor in cs
        nvertices = max(nvertices, length(neighbor.vertices))
    end
    return StaticParabolicNgon(c, p, max(8, 2nvertices + 2))
end

"""
    pmof(initial, fraction, first_moment, curvature, cell; origin, workspace)

Parabolic MOF (PMOF) reconstruction.  As in the restricted `Q²_κ` method in
the paper, the curvature is an explicit input; volume fraction and first
moment alone do not determine its sign.  `origin` only sets the local
coordinate frame and defaults to `centroid(cell)`.
"""
function pmof(p0::PlanarHS{2}, α::Real, first_moment, curvature::Quantity, c::Ngon;
    origin::Point=centroid(c),
    workspace::Union{Nothing,StaticParabolicNgon}=nothing,
    xatol::Real=√eps(Float64),
    maxiters::Integer=25,
    step_max::Real=.5,
    verbose::Bool=false,
)
    _validate_fraction(α)
    cell_area = abs(smeasure(c))
    target_area = α * cell_area
    moment_scale = cell_area * sqrt(cell_area)
    θ0 = normal_to_angle(p0.𝛈)
    seed = Parabola(angle_to_normal(θ0), zero(c.vertices[1].coords.x - origin.coords.x), curvature, origin)
    out = isnothing(workspace) ? _parabolic_workspace(c, seed) : workspace
    cost_and_derivative(θ) = begin
        p = _parabola_with_area(θ, curvature, target_area, c, origin; workspace=out)
        intersect!(out, c, p)
        _, reconstructed_moment = moments(out)
        difference = (reconstructed_moment - first_moment) / moment_scale
        dshift_angle, _ = _parabolic_shift_derivatives(out)
        derivative = _parabolic_first_moment_angle_derivative(out, dshift_angle) / moment_scale
        return sum(abs2, difference), 2 * dot(difference, derivative)
    end
    θ = brent_min(cost_and_derivative, θ0;
        xatol=xatol, maxiters=maxiters, step_max=step_max, verbose=verbose)
    return _parabola_with_area(θ, curvature, target_area, c, origin; workspace=out)
end

function _parabolic_lvira_costfun(p::Parabola, fractions, cells, areas, central, workspace)
    intersect!(workspace, central, p)
    dshift_angle, dshift_curvature = _parabolic_shift_derivatives(workspace)
    cost = zero(smeasure(workspace) / areas[1])
    dangle = zero(cost)
    dcurvature = zero(cost * central.vertices[1].coords.x)
    for (α, c, area) in zip(fractions, cells, areas)
        intersect!(workspace, c, p)
        error = smeasure(workspace) / area - α
        weight = 1 / (α * (1 - α) + 1e-2)
        cost += weight * error^2
        dvolume_angle, dvolume_curvature = _parabolic_volume_derivatives(workspace,
            dshift_angle, dshift_curvature)
        dangle += 2 * weight * error * dvolume_angle / area
        dcurvature += 2 * weight * error * dvolume_curvature / area
    end
    return cost, dangle, dcurvature
end

"""
    plvira(initial, fraction, curvature, cell, neighbor_fractions, neighbor_cells;
           cmeasures=smeasure.(neighbor_cells), origin, workspace)

Parabolic LVIRA (PLVIRA) reconstruction with a supplied curvature.  It fits
the volume fractions in the neighboring stencil while enforcing the central
cell fraction exactly.
"""
function plvira(p0::PlanarHS{2}, α::Real, curvature::Quantity, c::Ngon,
    αs::AbstractArray{<:Real}, cs::Union{SubDomain,AbstractArray{<:Ngon}};
    cmeasures::AbstractArray{<:Quantity}=smeasure.(cs),
    origin::Point=centroid(c),
    workspace::Union{Nothing,StaticParabolicNgon}=nothing,
    xatol::Real=√eps(Float64),
    maxiters::Integer=25,
    step_max::Real=.5,
    verbose::Bool=false,
)
    _validate_fraction(α)
    length(αs) == length(cs) == length(cmeasures) || throw(DimensionMismatch(
        "neighbor fractions, cells, and cell areas must have the same length"))
    target_area = α * smeasure(c)
    θ0 = normal_to_angle(p0.𝛈)
    seed = Parabola(angle_to_normal(θ0), zero(c.vertices[1].coords.x - origin.coords.x), curvature, origin)
    out = isnothing(workspace) ? _parabolic_workspace(c, seed, cs) : workspace
    cost_and_derivative(θ) = begin
        p = _parabola_with_area(θ, curvature, target_area, c, origin; workspace=out)
        cost, derivative, _ = _parabolic_lvira_costfun(p, αs, cs, cmeasures, c, out)
        return cost, derivative
    end
    θ = brent_min(cost_and_derivative, θ0;
        xatol=xatol, maxiters=maxiters, step_max=step_max, verbose=verbose)
    return _parabola_with_area(θ, curvature, target_area, c, origin; workspace=out)
end

"""
    prost(initial, fraction, cell, neighbor_fractions, neighbor_cells;
          curvature_bounds, cmeasures=smeasure.(neighbor_cells), origin)

Full quadratic (`Q²`) parabolic reconstruction.  PROST jointly searches the
normal and curvature that minimize the LVIRA volume mismatch.  Supplying
physical `curvature_bounds` is recommended when a stencil has more than one
plausible parabolic fit. The coupled solve is performed by `Optim.jl`'s
bounded L-BFGS implementation using the analytic cost gradient.
"""
function prost(p0::PlanarHS{2}, α::Real, c::Ngon,
    αs::AbstractArray{<:Real}, cs::Union{SubDomain,AbstractArray{<:Ngon}};
    cmeasures::AbstractArray{<:Quantity}=smeasure.(cs),
    curvature_bounds=nothing,
    origin::Point=centroid(c),
    workspace::Union{Nothing,StaticParabolicNgon}=nothing,
    maxiters::Integer=100,
)
    _validate_fraction(α)
    length(αs) == length(cs) == length(cmeasures) || throw(DimensionMismatch(
        "neighbor fractions, cells, and cell areas must have the same length"))
    target_area = α * smeasure(c)
    scale = sqrt(abs(smeasure(c)))
    raw_κbounds = isnothing(curvature_bounds) ? (-8 / scale, 8 / scale) : curvature_bounds
    κbounds = (float(raw_κbounds[1]), float(raw_κbounds[2]))
    κbounds[1] < κbounds[2] || throw(ArgumentError("curvature_bounds must be increasing"))
    θ0, κ0 = normal_to_angle(p0.𝛈), zero(κbounds[1])
    seed = Parabola(angle_to_normal(θ0), zero(c.vertices[1].coords.x - origin.coords.x), κ0, origin)
    out = isnothing(workspace) ? _parabolic_workspace(c, seed, cs) : workspace
    function cost_and_gradient(x)
        θ, scaled_curvature = x
        κ = scaled_curvature / scale
        p = _parabola_with_area(θ, κ, target_area, c, origin; workspace=out)
        cost, dangle, dcurvature = _parabolic_lvira_costfun(p, αs, cs, cmeasures, c, out)
        return cost, dangle, dcurvature / scale
    end
    objective(x) = first(cost_and_gradient(x))
    function gradient!(storage, x)
        _, dangle, dcurvature = cost_and_gradient(x)
        storage[1] = dangle
        storage[2] = dcurvature
        return storage
    end
    lower = [θ0 - π, ustrip(κbounds[1] * scale)]
    upper = [θ0 + π, ustrip(κbounds[2] * scale)]
    initial = [θ0, clamp(0.0, lower[2], upper[2])]
    result = Optim.optimize(objective, gradient!, lower, upper, initial,
        Optim.Fminbox(Optim.LBFGS(linesearch=Optim.LineSearches.BackTracking())),
        Optim.Options(iterations=maxiters))
    θ, scaled_curvature = Optim.minimizer(result)
    return _parabola_with_area(θ, scaled_curvature / scale, target_area, c, origin; workspace=out)
end
