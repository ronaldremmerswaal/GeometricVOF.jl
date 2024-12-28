"""
    V = volume(Φ::Function, c::Ngon)

Compute the volume V occupied by the reference phase. I.e. compute the volume of the set
    𝕊 = {𝐱 ∈ c | Φ(𝐱) ≤ 0}.
"""
function volume(Φ::Function, c::Ngon)
    verts = c.vertices

    Φp(p::Point) = Φ(to(p)...)
    T = typeof(Φp(verts[1]))

    if T == Bool
        # Here Φ(𝐱) == true ⟺ 𝐱 ∈ 𝕊, and we must use bisection on a shifted function
        root_finder = Roots.Bisection()
        Φrf(p) = Φp(p) - .5
    else
        root_finder = Roots.Brent()
        Φrf = Φp
    end

    Φverts = Φrf.(verts)
    inside_verts = Φverts ≤ 0

    # Quick return for trivial cases
    if all(inside_verts)
        return measure(c)
    elseif all(!inside_verts)
        return 0u"m^2"
    end

    # TODO: compute polygonal approx of reference phase inside c
end
