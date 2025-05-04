tangent(s::Segment{𝔼{2}}) = tangent!(zeros(2), s)

function tangent!(𝛕::AbstractVector, s::Segment{𝔼{2}})
    𝛕[1] = ustrip(s.vertices[2].coords.x - s.vertices[1].coords.x)
    𝛕[2] = ustrip(s.vertices[2].coords.y - s.vertices[1].coords.y)
    normalize!(𝛕)
end


normal_x(s::Segment{𝔼{2}}) = ustrip(s.vertices[2].coords.y - s.vertices[1].coords.y) / √(ustrip(s.vertices[2].coords.x - s.vertices[1].coords.x)^2 + ustrip(s.vertices[2].coords.y - s.vertices[1].coords.y)^2)
normal_y(s::Segment{𝔼{2}}) = -ustrip(s.vertices[2].coords.x - s.vertices[1].coords.x) / √(ustrip(s.vertices[2].coords.x - s.vertices[1].coords.x)^2 + ustrip(s.vertices[2].coords.y - s.vertices[1].coords.y)^2)

normal(s::Segment{𝔼{2}}) = normal!(zeros(2), s)
function normal!(𝛈::AbstractVector, s::Segment{𝔼{2}})
    tangent!(𝛈, s)
    𝛈[1], 𝛈[2] = 𝛈[2], -𝛈[1]
    𝛈
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

function donating_region(s::Segment{𝔼{2}}, velo::AbstractVector,
    dt::Quantity; α::Union{Nothing, Quantity}=nothing)

    out = StaticNgon(eltype(s.vertices))
    Ngon(donating_region!(out, s, velo, dt; α=α))
end

function donating_region!(out::StaticNgon, s::Segment{𝔼{2}}, velo::AbstractVector,
    dt::Quantity; α::Union{Nothing, Quantity}=nothing)

    v1, v2 = vertices(s)

    v1_pre = v1 - Vec(dt * velo[1][1], dt * velo[1][2])
    v2_pre = v2 - Vec(dt * velo[2][1], dt * velo[2][2])

    out.vertices[1] = v1
    out.vertices[2] = v2
    out.vertices[3] = v2_pre
    out.vertices[4] = v1_pre
    out.nr_verts = 4

    if !isnothing(α)
        s_pre = Segment(v1_pre, v2_pre)
        n_x = normal_x(s_pre)
        n_y = normal_y(s_pre)

        volume_err = smeasure(out) - α
        scaling = 2 *  volume_err / Meshes.measure(s_pre)

        vmid_pre_x = (v1_pre.coords.x + v2_pre.coords.x) / 2 + n_x * scaling
        vmid_pre_y = (v1_pre.coords.y + v2_pre.coords.y) / 2 + n_y * scaling

        out.vertices[4] = Point(vmid_pre_x, vmid_pre_y)
        out.vertices[5] = v1_pre
        out.nr_verts = 5
    end

    return out
end
