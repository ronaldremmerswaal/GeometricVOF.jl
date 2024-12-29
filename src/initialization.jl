import Meshes.measure

"""
    M = measure(Φ::Function, c::Ngon)

Compute the measure M occupied by the reference phase. I.e. compute the measure of the set
    𝕊 = {𝐱 ∈ c | Φ(𝐱) ≤ 0},
or
    𝕊 = {𝐱 ∈ c | Φ(𝐱)}
if Φ is boolean.
"""
function Meshes.measure(Φ::Function, c::Ngon)
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
    𝛈 = Vec(𝛕[2], -𝛕[1]) # Outward pointing normal (in hf direction)

    # The levelset function in local coordinates on τ ∈ [0, 1]
    Φl(τ, η) = Φ(v2 + τ * 𝛕 + η * 𝛈)

    # For each τ, the value of the height-function results from solving a rootfinding problem
    hf(τ) = find_zero(η -> Φl(τ, η), (-ηmax, ηmax), Method())

    τ_gl, weight_gl = gausslegendre(16) # TODO: precompute?
    τ_gl = (τ_gl .+ 1)/2
    h = norm(𝛕)

    return -h^2 * (hf.(τ_gl) ⋅ weight_gl) / 2
end
