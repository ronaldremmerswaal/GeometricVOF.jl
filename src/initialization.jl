
"""
    M = smeasure(Φ::Function, c::Ngon)

Compute the measure M occupied by the reference phase. I.e. compute the measure of the set
    𝕊 = {𝐱 ∈ c | Φ(𝐱) ≤ 0},
or
    𝕊 = {𝐱 ∈ c | Φ(𝐱)}
if Φ is boolean.
In rare cases an approximation error might occur, which can be approximated as (relative to
the diameter H of the Ngon):
    H κ / (8 * nref^2)
where κ is a bound on the curvature of the interface and h = H / nref. Here, we assume
κ < 1 / h.
"""
function smeasure(Φ::F, c::Ngon; nref::Int=8, workspaces=(StaticNgon(c, nref*8), StaticNgon(c, nref*8)), gl_data=gausslegendre(16)) where {F}
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

    ref_c, workspace = workspaces

    refine_edges!(ref_c, c, nref=nref)  # Optional step that allows for detection of multiple
                                    # intersections per edge

    Φverts = Φrf.(view(ref_c.vertices, 1:ref_c.nr_verts))
    inside_verts = Φverts .≤ zero(T)

    # Quick return for trivial cases
    if all(inside_verts)
        return smeasure(ref_c)
    elseif all(.!inside_verts)
        return 0u"m^2"
    end

    workspace, edge_is_hf = ngon_approx!(workspace, Φrf, ref_c, inside_verts, Method)
    M = smeasure(workspace)
    for edx = 1 : workspace.nr_verts
        if edge_is_hf[edx]
            # We correct the measure M for each edge which is actually a height-function
            v1 = workspace.vertices[edx]
            v2 = workspace.vertices[mod1(edx + 1, workspace.nr_verts)]

            M += hf_measure(Φrf, v1, v2, Method, gl_data)
        end
    end

    return M
end

function ngon_approx!(Φc_approx::StaticNgon{N}, Φ::F, c::StaticNgon, inside_verts::AbstractVector{Bool}, Method) where {N, F}
    verts = c.vertices

    vdx = findfirst(inside_verts)

    edge_is_hf = MVector{N, Bool}(undef)

    new_idx = 0
    for _ = 1 : c.nr_verts
        ndx = mod1(vdx + 1, c.nr_verts)

        if inside_verts[vdx]
            new_idx += 1
            Φc_approx.vertices[new_idx] = verts[vdx]
            edge_is_hf[new_idx] = false
        end

        if inside_verts[vdx] != inside_verts[ndx]
            # Compute the new vertex
            v1 = verts[vdx]
            v2 = verts[ndx]
            v(s) = v1 + (v2 - v1) * s

            new_idx += 1
            Φc_approx.vertices[new_idx] = v(find_zero(Φ ∘ v, (0, 1), Method()))
            edge_is_hf[new_idx] = inside_verts[vdx]
        end

        vdx = ndx
    end

    Φc_approx.nr_verts = new_idx

    return Φc_approx, edge_is_hf
end

function hf_measure(Φ::F, v1::Point, v2::Point, Method, gl_data) where {F}
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

    τ_gl, weight_gl = gl_data

    ans = 0
    for (τ, weight) ∈ zip(τ_gl, weight_gl)
        ans += hf((τ + 1)/2) * weight
    end

    return -h^2 * ans / 2
end

function refine_edges!(ref_c::StaticNgon, c::Ngon; nref::Int=2)
    @assert nref ≥ 1 "Refinement nref (given by $nref) must be a positive integer"

    vs = vertices(c)

    new_idx = 0
    for vdx ∈ eachindex(vs)
        ndx = mod1(vdx + 1, length(vs))

        v = vs[vdx]
        dir = vs[ndx] - v

        new_idx += 1
        ref_c.vertices[new_idx] = v
        for rdx = 1 : nref - 1
            new_idx += 1
            ref_c.vertices[new_idx] = v + (rdx / nref) * dir
        end
    end

    ref_c.nr_verts = new_idx

    return ref_c
end
