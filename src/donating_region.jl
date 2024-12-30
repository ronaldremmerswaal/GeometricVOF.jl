tangent(s::Segment{𝔼{2}}) = normalize(s.vertices[2] - s.vertices[1])

function normal(s::Segment{𝔼{2}})
    𝐭 = tangent(s)
    [-𝐭[2], 𝐭[1]]
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
    if isnothing(α)
        # No volume constraint
        return poly
    end

    s_pre = Segment(v1_pre, v2_pre)
    n = normal(s_pre)
    volume_err = α - measure(poly)
    dn = n * (2 *  volume_err / measure(s_pre))
    vmid_pre = v1_pre + (v2_pre - v1_pre) / 2 + Vec(dn...)

    Pentagon(v1, v2, v2_pre, vmid_pre, v1_pre)
end
