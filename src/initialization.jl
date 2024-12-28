"""
    V = volume(Φ::Function, c::Ngon)

Compute the volume V occupied by the reference phase. I.e. compute the volume of the set
    𝕊 = {𝐱 ∈ c | Φ(𝐱) ≤ 0},
or
    𝕊 = {𝐱 ∈ c | Φ(𝐱)}
if Φ is boolean.
"""
function volume(Φ::Function, c::Ngon)
    verts = c.vertices

    Φp(p::Point) = Φ(to(p)...)
    T = typeof(Φp(verts[1]))

    if T == Bool
        # Here Φ(𝐱) == true ⟺ 𝐱 ∈ 𝕊, and we must use bisection on a shifted function
        Method = Roots.Bisection
        Φrf(p) = Φp(p) - .5
    else
        Method = Roots.Brent
        Φrf = Φp
    end

    Φverts = Φrf.(verts)
    inside_verts = Φverts .≤ zero(T)

    # Quick return for trivial cases
    if all(inside_verts)
        return measure(c)
    elseif all(.!inside_verts)
        return 0u"m^2"
    end

    # TODO: compute polygonal approx of reference phase inside c
    ϕc_approx, edge_is_hf = ngon_approx(Φrf, c, inside_verts, Method)
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

    return new_verts, edge_is_hf
end
