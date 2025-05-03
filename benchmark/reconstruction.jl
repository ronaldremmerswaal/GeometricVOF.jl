using GeometricVOF

function reconstruction_prep(;N=2^6)
    R = 0.25u"m"
    Φ(x, y) = (x - .5u"m")^2 + (y - .5u"m")^2 - R^2 * (1 + .1cos(.1 + 5atan(y - .5u"m", x - .5u"m")))^2

    hx = 1.333333u"m" / N
    hy = 0.966666u"m" / N
    mesh = CartesianGrid((N, N), (0., 0.), (hx, hy))

    inds = LinearIndices(size(mesh))
    αs = reshape([measure(Φ, c) / measure(c) for c ∈ mesh], (N, N))

    return mesh, inds, αs
end

function reconstruction_benchmark(mesh, inds, αs; α_tol=1E-8)
    recons = PlanarHS[]
    for i = 2 : size(inds, 1)-1, j = 2 : size(inds, 2)-1
        if αs[i, j] < α_tol || αs[i, j] > 1 - α_tol
            continue
        end
        c = mesh[i, j]
        xc = centroid(c)
        θ0 = atan(xc.coords.y, xc.coords.x) + .1
        p0 = PlanarHS(GeometricVOF.angle_to_normal(θ0), 0u"m")

        p_recon = reconstruct(p0, αs[i, j], c, view(mesh, inds[i-1:i+1, j-1:j+1][:]), view(αs, i-1:i+1, j-1:j+1))
        push!(recons, p_recon)
    end

    recons
end
