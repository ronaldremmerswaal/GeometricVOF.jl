import Meshes.measure

"""
    M = measure(Φ::Function, c::Ngon)

Compute the measure M occupied by the reference phase. I.e. compute the measure of the set
    𝕊 = {𝐱 ∈ c | Φ(𝐱) ≤ 0},
or
    𝕊 = {𝐱 ∈ c | Φ(𝐱)}
if Φ is boolean.
In rare cases an approximation error might occur, which can be approximated as (relative to
the diameter H of the Ngon):
    π / (2 κ^2)             if κ > 1 / h
    H κ / (8 * nref^2)      if κ < 1 / h
where κ is a bound on the curvature of the interface and h = H / nref
"""
function Meshes.measure(Φ::Function, c::Ngon; nref::Int=8)
    Φp(p::Point) = Φ(p.coords.x, p.coords.y)
    T = typeof(Φp(c.vertices[1]))

    if T == Bool
        # Here Φ(𝐱) == true ⟺ 𝐱 ∈ 𝕊, and we must use bisection on a shifted function
        Method = Roots.Bisection
        Φrf(p) = .5 - Φp(p)
    else
        Method = Roots.Brent
        Φrf = Φp
    end

    c = refine_edges(c, nref=8)   # Optional step that allows for detection of multiple
                                    # intersections per edge

    Φverts = Φrf.(c.vertices)
    inside_verts = Φverts .≤ zero(T)

    # Quick return for trivial cases
    if all(inside_verts)
        return measure(c)
    elseif all(.!inside_verts)
        return 0u"m^2"
    end

    Φc_approx, edge_is_hf = ngon_approx(Φrf, c, inside_verts, Method)
    M = measure(Φc_approx)
    for (edx, is_hf) = enumerate(edge_is_hf)
        if is_hf
            # We correct the measure M for each edge which is actually a height-function
            v1 = Φc_approx.vertices[edx]
            v2 = Φc_approx.vertices[edx == length(edge_is_hf) ? 1 : edx + 1]

            M += hf_measure(Φrf, v1, v2, Method)
        end
    end

    return M
end

function ngon_approx(Φ::Function, c::Ngon, inside_verts::AbstractVector{Bool}, Method)
    verts = c.vertices

    first_vdx = findfirst(inside_verts)
    vdx = first_vdx

    new_verts = Vector{eltype(verts)}()
    edge_is_hf = Vector{Bool}()

    for _ = 1 : length(verts)
        next_vdx = vdx == length(verts) ? 1 : vdx + 1

        if inside_verts[vdx]
            push!(new_verts, verts[vdx])
            push!(edge_is_hf, false)
        end

        if inside_verts[vdx] != inside_verts[next_vdx]
            # Compute the new vertex
            v1 = verts[vdx]
            v2 = verts[next_vdx]
            v(s) = v1 + (v2 - v1) * s

            push!(new_verts, v(find_zero(Φ ∘ v, (0, 1), Method())))
            push!(edge_is_hf, inside_verts[vdx])
        end

        vdx = next_vdx
    end

    return Ngon(new_verts...), edge_is_hf
end

function hf_measure(Φ::Function, v1::Point, v2::Point, Method)
    ηmax = 1    # Length-scale (relative to norm(v1 - v2)) used in find_zero

    𝛕 = v1 - v2
    h = norm(𝛕)
    if h < 10eps()u"m"
        return 0u"m^2"
    end
    𝛈 = Vec(𝛕[2], -𝛕[1]) # Outward pointing normal (in hf direction)

    # The levelset function in local coordinates on τ ∈ [0, 1]
    Φl(τ, η) = Φ(v2 + τ * 𝛕 + η * 𝛈)

    # For each τ, the value of the height-function results from solving a rootfinding problem
    hf(τ) = find_zero(η -> Φl(τ, η), (-ηmax, ηmax), Method())

    τ_gl, weight_gl = gausslegendre(16) # TODO: precompute?
    τ_gl = (τ_gl .+ 1)/2

    return -h^2 * (hf.(τ_gl) ⋅ weight_gl) / 2
end

function refine_edges(c::Ngon; nref::Int=2)
    @assert nref ≥ 1 "Refinement nref (given by $nref) must be a positive integer"
    if nref == 1
        return c
    end
    vs = vertices(c)

    rvs = Vector{eltype(vs)}()
    for vdx ∈ eachindex(vs)
        ndx = vdx == length(vs) ? 1 : vdx + 1

        v = vs[vdx]
        dir = vs[ndx] - vs[vdx]
        push!(rvs, v)
        for rdx = 1 : nref - 1
            push!(rvs, v + (rdx / nref) * dir)
        end
    end

    return Ngon(rvs...)
end
