"""
    reconstruct(initial, fraction, cell, neighbor_fractions, neighbor_cells;
                cell_areas=smeasure.(neighbor_cells), workspace, shift_workspace)

Reconstruct a planar interface in `cell` with LVIRA. `fraction` must be
strictly between zero and one. The allocating method is intended for ordinary
use; reuse the optional workspaces, or use `reconstruct!`, in a cell loop.
"""
function reconstruct(p0::PlanarHS{2}, α_central::T, c_central::Ngon,
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

"""
    reconstruct!(out, initial, fraction, cell, neighbor_fractions, neighbor_cells;
                 cell_areas=smeasure.(neighbor_cells), workspace, shift_workspace)

Allocation-free LVIRA reconstruction. `out`, `cell`, `neighbor_cells`, and the
workspace must be fixed-capacity `StaticNgon`s with enough room for clipped
intermediate polygons. Returns `out`.
"""
function reconstruct!(out::StaticNgon, p0::PlanarHS{2}, α_central::T,
    c_central::StaticNgon{N, P}, αs::AbstractVector{T}, cs::AbstractVector{<:StaticNgon},
    cmeasures::AbstractVector{Q}=smeasure.(cs);
    verbose::Bool=false,
    xatol::Real=√eps(T),
    workspace::StaticNgon=StaticNgon(P, N + 2),
    shift_workspace::AbstractVector{<:Real}=MVector{N + 2, Float64}(undef),
) where {T<:Real, Q<:Quantity, N, P<:Point}
    _validate_reconstruction_inputs(α_central, αs, cs, cmeasures)
    ref_vol = smeasure(c_central) * α_central

    wrapped_costfun(θ::Real) = lvira_costfun(PlanarHS(θ, ref_vol, c_central; workspace=workspace, shift_workspace=shift_workspace), cs, αs, cmeasures, c_central, workspace=workspace)

    θ0 = GeometricVOF.normal_to_angle(p0.𝛈)
    θ = brent_min(wrapped_costfun, θ0; xatol=xatol, maxiters=25, step_max=.5, verbose=verbose)

    p = PlanarHS(θ, ref_vol, c_central; workspace=workspace, shift_workspace=shift_workspace)
    intersect!(out, c_central, p)
    return out
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

@inline _moment_cost(reference, reconstructed) = sum(abs2, reference - reconstructed)

"""Minimize a scalar periodic objective without relying on a derivative."""
function _periodic_minimum(f, θ0::Real; samples::Integer=32, iterations::Integer=36)
    samples ≥ 3 || throw(ArgumentError("at least three angular samples are required"))
    Δθ = 2π / samples
    θbest = θ0 - π
    fbest = f(θbest)
    for index in 1:samples-1
        θ = θ0 - π + index * Δθ
        value = f(θ)
        if value < fbest
            θbest, fbest = θ, value
        end
    end
    return _golden_minimum(f, θbest - Δθ, θbest + Δθ; iterations=iterations)
end

function _golden_minimum(f, left, right; iterations::Integer=36)
    ratio = (sqrt(5) - 1) / 2
    x1 = right - ratio * (right - left)
    x2 = left + ratio * (right - left)
    f1, f2 = f(x1), f(x2)
    for _ in 1:iterations
        if f1 ≤ f2
            right, x2, f2 = x2, x1, f1
            x1 = right - ratio * (right - left)
            f1 = f(x1)
        else
            left, x1, f1 = x1, x2, f2
            x2 = left + ratio * (right - left)
            f2 = f(x2)
        end
    end
    return (left + right) / 2
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
    samples::Integer=32,
    iterations::Integer=36,
)
    _validate_fraction(α)
    target_area = α * smeasure(c)
    cost(θ) = begin
        p = PlanarHS(θ, target_area, c; workspace=workspace, shift_workspace=shift_workspace)
        _, reconstructed_moment = moments(p, c; workspace=workspace)
        _moment_cost(first_moment, reconstructed_moment)
    end
    θ = _periodic_minimum(cost, normal_to_angle(p0.𝛈); samples=samples, iterations=iterations)
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
    samples::Integer=32,
    iterations::Integer=36,
)
    _validate_fraction(α)
    target_area = α * smeasure(c)
    θ0 = normal_to_angle(p0.𝛈)
    seed = Parabola(angle_to_normal(θ0), zero(c.vertices[1].coords.x - origin.coords.x), curvature, origin)
    out = isnothing(workspace) ? _parabolic_workspace(c, seed) : workspace
    cost(θ) = begin
        p = _parabola_with_area(θ, curvature, target_area, c, origin; workspace=out)
        _, reconstructed_moment = moments(p, c; workspace=out)
        _moment_cost(first_moment, reconstructed_moment)
    end
    θ = _periodic_minimum(cost, θ0; samples=samples, iterations=iterations)
    return _parabola_with_area(θ, curvature, target_area, c, origin; workspace=out)
end

function _parabolic_lvira_cost(p::Parabola, fractions, cells, areas, workspace)
    cost = zero(smeasure(p, cells[1]; workspace=workspace) / areas[1])
    for (α, c, area) in zip(fractions, cells, areas)
        error = smeasure(p, c; workspace=workspace) / area - α
        weight = 1 / (α * (1 - α) + 1e-2)
        cost += weight * error^2
    end
    return cost
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
    samples::Integer=32,
    iterations::Integer=36,
)
    _validate_fraction(α)
    length(αs) == length(cs) == length(cmeasures) || throw(DimensionMismatch(
        "neighbor fractions, cells, and cell areas must have the same length"))
    target_area = α * smeasure(c)
    θ0 = normal_to_angle(p0.𝛈)
    seed = Parabola(angle_to_normal(θ0), zero(c.vertices[1].coords.x - origin.coords.x), curvature, origin)
    out = isnothing(workspace) ? _parabolic_workspace(c, seed, cs) : workspace
    cost(θ) = begin
        p = _parabola_with_area(θ, curvature, target_area, c, origin; workspace=out)
        _parabolic_lvira_cost(p, αs, cs, cmeasures, out)
    end
    θ = _periodic_minimum(cost, θ0; samples=samples, iterations=iterations)
    return _parabola_with_area(θ, curvature, target_area, c, origin; workspace=out)
end

"""
    prost(initial, fraction, cell, neighbor_fractions, neighbor_cells;
          curvature_bounds, cmeasures=smeasure.(neighbor_cells), origin)

Full quadratic (`Q²`) parabolic reconstruction.  PROST jointly searches the
normal and curvature that minimize the LVIRA volume mismatch.  Supplying
physical `curvature_bounds` is recommended when a stencil has more than one
plausible parabolic fit.
"""
function prost(p0::PlanarHS{2}, α::Real, c::Ngon,
    αs::AbstractArray{<:Real}, cs::Union{SubDomain,AbstractArray{<:Ngon}};
    cmeasures::AbstractArray{<:Quantity}=smeasure.(cs),
    curvature_bounds=nothing,
    origin::Point=centroid(c),
    workspace::Union{Nothing,StaticParabolicNgon}=nothing,
    angle_samples::Integer=24,
    curvature_samples::Integer=9,
    iterations::Integer=24,
)
    _validate_fraction(α)
    length(αs) == length(cs) == length(cmeasures) || throw(DimensionMismatch(
        "neighbor fractions, cells, and cell areas must have the same length"))
    target_area = α * smeasure(c)
    scale = sqrt(abs(smeasure(c)))
    raw_κbounds = isnothing(curvature_bounds) ? (-8 / scale, 8 / scale) : curvature_bounds
    κbounds = (float(raw_κbounds[1]), float(raw_κbounds[2]))
    κbounds[1] < κbounds[2] || throw(ArgumentError("curvature_bounds must be increasing"))
    θbest, κbest = normal_to_angle(p0.𝛈), zero(κbounds[1])
    seed = Parabola(angle_to_normal(θbest), zero(c.vertices[1].coords.x - origin.coords.x), κbest, origin)
    out = isnothing(workspace) ? _parabolic_workspace(c, seed, cs) : workspace
    function cost(θ, κ)
        p = _parabola_with_area(θ, κ, target_area, c, origin; workspace=out)
        return _parabolic_lvira_cost(p, αs, cs, cmeasures, out)
    end
    best_cost = cost(θbest, κbest)
    κgrid = range(κbounds[1], κbounds[2]; length=curvature_samples)
    for κ in κgrid
        θ = _periodic_minimum(θ -> cost(θ, κ), θbest;
            samples=angle_samples, iterations=iterations)
        value = cost(θ, κ)
        if value < best_cost
            θbest, κbest, best_cost = θ, κ, value
        end
    end
    # Coordinate refinement keeps the public API compact while providing a
    # dependable full-Q² solve without derivatives of the clipping algorithm.
    for _ in 1:3
        θbest = _periodic_minimum(θ -> cost(θ, κbest), θbest;
            samples=angle_samples, iterations=iterations)
        κbest = _golden_minimum(κ -> cost(θbest, κ), κbounds[1], κbounds[2];
            iterations=iterations)
    end
    return _parabola_with_area(θbest, κbest, target_area, c, origin; workspace=out)
end
