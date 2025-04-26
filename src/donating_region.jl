tangent(s::Segment{𝔼{2}}) = normalize(s.vertices[2] - s.vertices[1])

function normal(s::Segment{𝔼{2}})
    𝐭 = tangent(s)
    [-𝐭[2], 𝐭[1]]
end

"""
    smeasure(p)

Signed measure of Ngon `p`.
"""
smeasure(p::Ngon) = smeasure(p.vertices)

"""
    smeasure(p)

Signed measure of PolyArea `p`.
"""
function smeasure(p::PolyArea)
    M = 0u"m^2"
    for ring ∈ p.rings
        M += smeasure(ring.vertices)
    end
    M
end

function smeasure(verts::AbstractVector{<:Point{𝔼{2}}})
    M = 0u"m^2"
    for vdx ∈ eachindex(verts)
        v1 = verts[vdx]
        v2 = verts[vdx==length(verts) ? 1 : vdx + 1]
        M += v1.coords.x * v2.coords.y - v1.coords.y * v2.coords.x
    end
    M /= 2
end

function smeasure(p::PlanarHS{2}, c::Ngon)
    cp = c ∩ p
    isnothing(cp) ? 0u"m^2" : smeasure(cp)
end

function donating_region(s::Segment{𝔼{2}}, u::Function, dt::Quantity;
    α::Union{Nothing, Quantity}=nothing)
    velo = [u(to(s.vertices[1])...), u(to(s.vertices[2])...)]
    donating_region(s, velo, dt, α=α)
end

function donating_region(s::Segment{𝔼{2}}, velo::AbstractVector{Vector{T}},
    dt::Quantity; α::Union{Nothing, Quantity}=nothing) where T<:Quantity

    v1, v2 = vertices(s)

    v1_pre = v1 - Vec(dt * velo[1][1], dt * velo[1][2])
    v2_pre = v2 - Vec(dt * velo[2][1], dt * velo[2][2])

    poly = Quadrangle(v1, v2, v2_pre, v1_pre)

    s_pre = Segment(v1_pre, v2_pre)
    n = normal(s_pre)

    if !isnothing(α)
        volume_err = α - smeasure(poly)
        dn = n * (2 *  volume_err / measure(s_pre))
        vmid_pre = v1_pre + (v2_pre - v1_pre) / 2 + Vec(dn...)

        poly = Pentagon(v1, v2, v2_pre, vmid_pre, v1_pre)
    end

    # We split the DR into at most two parts: one positively oriented, one negatively
    hs = PlanarHS(n, n ⋅ to(v1))
    (poly ∩ complement(hs), poly ∩ hs)
end
