using GeometricVOF
using GLMakie
using Meshes
using Unitful

# The interface is a flower shape
R0 = 0.3u"m"
R(θ) = R0 * (1 + 0.1 * sin(0.1 + 5θ))
Φ(x, y) = x^2 + y^2 - R(atan(y, x))^2

N = 32
h = 1 / N
mesh = CartesianGrid((N, N), (-.5, -.5), (h, h))

fig = Figure()
ax = Axis(fig[1, 1], aspect=1)

# Show grid lines
viz!(ax, mesh, showsegments=true, color=:white, segmentcolor=:gray)

αs = reshape([measure(Φ, c) for c ∈ mesh], (N, N))
for i = 2 : N-1, j = 2 : N-1
    # On each interior cell the interface is reconstructed and plotted
    recon = LVIRA(collect(mesh[i-1:i+1, j-1:j+1]), αs[i-1:i+1, j-1:j+1])
    c = mesh[i, j]
    xc = centroid(c)
    θ0 = atan(xc.coords.y, xc.coords.x) + .1
    p0 = PlanarHS{2}(GeometricVOF.angle_to_normal(θ0), 0u"m")
    p_recon = reconstruct(recon, c, αs[i, j], p0)

    cp = c ∩ p_recon
    if !isnothing(cp)
        viz!(ax, cp, alpha=.8)
    end
end

# Plot the exact interface
θs = range(0, 2π, 100)
lines!(ax, ustrip(R.(θs) .* cos.(θs)), ustrip(R.(θs) .* sin.(θs)), linewidth=4)

fig
