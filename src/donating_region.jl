function donating_region(s::Segment, velo::AbstractVector{Vec{2}}, dt::Quantity;
    α::Union{Nothing, Quantity}=nothing)

    v1, v2 = vertices(s)

    v1_pre = v1 - dt * velo[1]
    v2_pre = v2 - dt * velo[2]

    poly = Quadrangle(v1, v2, v2_pre, v1_pre)
    if isnothing(α)
        # No volume constraint
        return poly
    end

    avec = area_vector(s) # TODO
    face_length_pre = norm(avec)
    n = normalize!(avec)
    volume_err = measure(poly) - α
    vmid_pre = (v1_pre + v2_pre) / 2 + n * (2 *  volume_err / face_length_pre)

    return Quadrangle(v1, v2, v2_pre, vmid_pre, v1_pre)
end
