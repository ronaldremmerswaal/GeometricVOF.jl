
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
function smeasure(Φ::Function, c::Ngon; nref::Int=8, workspace::StaticNgon=StaticNgon(c, nref*8), gl_data=gausslegendre(16))
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

    c = refine_edges(c, nref=nref)  # Optional step that allows for detection of multiple
                                    # intersections per edge

    Φverts = Φrf.(c.vertices)
    inside_verts = Φverts .≤ zero(T)

    # Quick return for trivial cases
    if all(inside_verts)
        return smeasure(c)
    elseif all(.!inside_verts)
        return 0u"m^2"
    end

    workspace, edge_is_hf = ngon_approx!(workspace, Φrf, c, inside_verts, Method)
    M = smeasure(workspace)
    for edx = 1 : workspace.nr_verts
        is_hf = edge_is_hf[edx]
        if is_hf
            # We correct the measure M for each edge which is actually a height-function
            v1 = workspace.vertices[edx]
            v2 = workspace.vertices[mod1(edx + 1, workspace.nr_verts)]

            M += hf_measure(Φrf, v1, v2, Method, gl_data)
        end
    end

    return M
end

function ngon_approx!(Φc_approx::StaticNgon{N}, Φ::Function, c::Ngon, inside_verts::AbstractVector{Bool}, Method) where N
    verts = c.vertices

    first_vdx = findfirst(inside_verts)
    vdx = first_vdx

    edge_is_hf = MVector{N, Bool}(undef)

    new_idx = 0
    for _ = 1 : length(verts)
        next_vdx = vdx == length(verts) ? 1 : vdx + 1

        if inside_verts[vdx]
            new_idx += 1
            Φc_approx.vertices[new_idx] = verts[vdx]
            edge_is_hf[new_idx] = false
        end

        if inside_verts[vdx] != inside_verts[next_vdx]
            # Compute the new vertex
            v1 = verts[vdx]
            v2 = verts[next_vdx]
            v(s) = v1 + (v2 - v1) * s

            new_idx += 1
            Φc_approx.vertices[new_idx] = v(find_zero(Φ ∘ v, (0, 1), Method()))
            edge_is_hf[new_idx] = inside_verts[vdx]
        end

        vdx = next_vdx
    end

    Φc_approx.nr_verts = new_idx

    return Φc_approx, edge_is_hf
end

function hf_measure(Φ::Function, v1::Point, v2::Point, Method, gl_data)
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
    τ_gl = (τ_gl .+ 1)/2

    ans = 0
    for (τ, weight) ∈ zip(τ_gl, weight_gl)
        ans += hf(τ) * weight
    end

    return -h^2 * ans / 2
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
