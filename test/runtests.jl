using GeometricVOF
using Meshes
using Unitful
using Test

@testset verbose=true "GeometricVOF.jl" begin
    @testset "initialize" begin
        # Trivial: all inside
        t = Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        for Φ ∈ [(x, y) -> x - 2u"m", (x, y) -> x ≤ 2u"m"]
            @test isapprox(smeasure(Φ, t), smeasure(t), rtol=10eps())
        end

        # Trivial: all outside
        for Φ ∈ [(x, y) -> x + 2u"m", (x, y) -> x ≤ -2u"m"]
            @test isapprox(smeasure(Φ, t), 0u"m^2", rtol=10eps())
        end

        # Linear
        for Φ ∈ [(x, y) -> x - 0.5u"m", (x, y) -> x ≤ 0.5u"m"]
            @test isapprox(smeasure(Φ, t), .5(1 - 1/4)u"m^2", rtol=10eps())
        end

        # Non-convex case
        nc = Ngon((0, 0), (1, 0), (.5, .5), (1, 1), (0, 1))
        for Φ ∈ [(x, y) -> x - .75u"m", (x, y) -> x ≤ .75u"m"]
            @test isapprox(smeasure(Φ, nc), (3/4 - 1/16)u"m^2", rtol=10eps())
        end

        # Monomial
        q = Quadrangle((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
        for p = 1 : 5
            for Φ ∈ [(x, y) -> y - (.5 + .25(x/u"m")^p)u"m", (x, y) -> y ≤ (.5 + .25(x/u"m")^p)u"m"]
                @test isapprox(smeasure(Φ, q), (1/2 + 1/(4(p+1)))u"m^2", rtol=10eps())
            end
        end

        # Circle
        ref(t) = refine(t, TriSubdivision())
        mesh = ref(ref(ref(q)))
        R = .25u"m"
        Φcircle(x, y) = (x - .5123u"m")^2 + (y - .45u"m")^2 - R^2
        M = 0.0u"m^2"
        for c ∈ mesh
            M += smeasure(Φcircle, c)
        end

        @test isapprox(M, π * R^2, rtol=10eps())

        # Edge case: colinear to edge, but all inside
        for Φl ∈ [(x, y) -> y - (1u"m" - x), (x, y) -> y ≤ 1u"m" - x]
            @test isapprox(smeasure(Φl, t), smeasure(t), rtol=10eps())
        end

        # Edge case: colinear to edge, all outside
        for Φl ∈ [(x, y) -> -y + (1u"m" - x), (x, y) -> y > 1u"m" - x]
            @test smeasure(Φl, t) < eps()u"m^2"
        end
    end

    @testset verbose=true "HalfSpace" begin
        @testset "intersect" begin
            t = Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))

            # Trivial: all inside
            tp = t ∩ PlanarHS([1, 0], 2u"m")
            @test tp == t

            # Trivial: all outside
            tp = t ∩ PlanarHS([-1, 0], -2u"m")
            @test isnothing(tp)

            # Quad to Quad
            c = Quadrangle((0, 0), (1, 0), (1, 1), (0, 1))
            p = PlanarHS([0, 1], .5u"m")
            cp = c ∩ p
            @test cp == Quadrangle((0, 0), (1, 0), (1, .5), (0, .5))

            # The third edge of the new polygon coincides with the HalfSpace
            e1 = Segment(cp.vertices[3], cp.vertices[4])
            @test GeometricVOF.normal(e1) == p.𝛈

            # Triangle to Quadrangle
            p = PlanarHS([1, 0], .5u"m")
            tp = t ∩ p
            @test tp == Quadrangle((0, 0), (.5, 0), (.5, .5), (0, 1))

            # The second edge of the new polygon coincides with the HalfSpace
            e1 = Segment(tp.vertices[2], tp.vertices[3])
            @test GeometricVOF.normal(e1) == p.𝛈

            # Non-convex case
            nc = Ngon((0, 0), (1, 0), (.5, .5), (1, 1), (0, 1))
            p = PlanarHS([1, 0], .75u"m")
            ncp = nc ∩ p
            @test ncp == Ngon((0, 0), (.75, 0), (.75, .25), (.5, .5), (.75, .75), (.75, 1), (0, 1))

            # The second edge of the new polygon coincides with the HalfSpace
            e1 = Segment(ncp.vertices[2], ncp.vertices[3])
            @test GeometricVOF.normal(e1) == p.𝛈

            # Edge case: colinear to edge, but all inside
            p = PlanarHS([1, 1] / √2, √.5u"m")
            tp = t ∩ p
            @test tp == t

            # Edge case: colinear to edge, but all outside
            tp = t ∩ complement(p)
            @test isnothing(tp)
        end

        @testset "shift" begin
            c = Triangle((0., 0.), (1., 0.), (0., 1.))
            @test shift(c, [1.0, 0.0], 0.21875u"m^2") == 0.25u"m"

            # Test if shift is the inverse of measure
            N = 7
            c = Ngon([(θ == 0 ? 0 : cos(θ), sin(θ)) for θ ∈ 2π*(0:N-1)/N]...) # Pacman

            M = 10
            for θ ∈ 2π*(0:M-1)/M
                𝛈 = GeometricVOF.SVector{2}(cos(θ), sin(θ))
                shift_min, shift_max = GeometricVOF.shift_extrema(c, 𝛈)
                for s_ref ∈ range(1.2shift_min, 1.2shift_max, length=M)
                    p = PlanarHS(𝛈, s_ref)

                    α = smeasure(p, c)
                    s_comp = shift(c, 𝛈, α)
                    if α == 0u"m^2"
                        @test s_comp ≤ shift_min
                    elseif α == smeasure(c)
                        @test s_comp ≥ shift_max || s_comp ≈ shift_max
                    else
                        @test isapprox(s_comp, s_ref, rtol=100eps())
                    end
                end
            end
        end
    end

    @testset "symmetric_difference (Δ)" begin
        # Exact: Δ = zero
        c = Quadrangle((0, 0), (1, 0), (1, 1), (0, 1))
        p = PlanarHS([0, 1], .5u"m")
        Φ(x, y) = y - .5u"m" # TODO also for Bool
        @test symmetric_difference(Φ, p, c) == 0u"m^2"

        # Simple: Δ is difference between planes
        p = PlanarHS([0, 1], .4u"m")
        @test isapprox(symmetric_difference(Φ, p, c), .1u"m^2", rtol=10eps())

        p = PlanarHS([0, -1], -.5u"m")
        @test isapprox(symmetric_difference(Φ, p, c), 1u"m^2", rtol=10eps())

        # Nontrivial: planar approximation to parabola (correct volume)
        f(x, ε) = .5u"m" + ε * (x - .5u"m")^2/u"m" - ε * u"m" / 12
        p = PlanarHS([0, 1], .5u"m")
        for ε ∈ 10. .^(-4 : -1)
            Φε(x, y) = y - f(x, ε)
            xs = 1/2 - 1/√12
            @test isapprox(symmetric_difference(Φε, p, c), ε * (4/3 * xs^3 - 2xs^2 +
                2/3 * xs - 1/24 + 1/8 - 1/12)u"m^2", rtol=10eps()/ε)
        end
    end

    @testset "reconstruction" begin
        # Test if linear interface is reconstructed exactly
        mesh = CartesianGrid((3, 3), (-.5, -.5), (1/3, 1/3))

        M = 10
        for θ ∈ 2π*(0:M-1)/M
            𝛈 = GeometricVOF.SVector{2}(cos(θ), sin(θ))
            shift_min, shift_max = GeometricVOF.shift_extrema(mesh[5], 𝛈)
            for shift_ref ∈ range(.9shift_min, .9shift_max, length=M)
                p_ref = PlanarHS(𝛈, shift_ref)

                # Initialize reference volumes
                αs = [smeasure(p_ref, c) / smeasure(c) for c ∈ mesh]
                p0 = PlanarHS(GeometricVOF.angle_to_normal(θ + 0.7), 0u"m")
                p_recon = reconstruct(p0, αs[5], mesh[5], αs, view(mesh, 1:9))

                @test isapprox(p_recon.shift, p_ref.shift, rtol=100eps())
                @test isapprox(p_recon.𝛈, p_ref.𝛈, rtol=100eps())
            end
        end

        # Test accuracy of reconstruction of flower shape
        R = 0.3u"m"
        Φ(x, y) = x^2 + y^2 - R^2 * (1 + .1sin(.1 + 5atan(y, x)))^2
        sd_errs = Vector{Float64}()u"m^2"
        for N ∈ 2 .^(4 : 6)
            h = 1 / N
            mesh = CartesianGrid((N, N), (-.5, -.5), (h, h))

            inds = LinearIndices(size(mesh))

            αs = reshape([smeasure(Φ, c) / smeasure(c) for c ∈ mesh], (N, N))
            sd_err = 0u"m^2"
            for i = 2 : N-1, j = 2 : N-1
                if αs[i, j] == 0 || αs[i, j] == 1
                    continue
                end
                c = mesh[i, j]
                xc = centroid(c)
                θ0 = atan(xc.coords.y, xc.coords.x) + .1
                p0 = PlanarHS(GeometricVOF.angle_to_normal(θ0), 0u"m")
                p_recon = reconstruct(p0, αs[i, j], c, view(αs, i-1:i+1, j-1:j+1), view(mesh, inds[i-1:i+1, j-1:j+1][:]))

                sd_err += symmetric_difference(Φ, p_recon, c)
            end

            push!(sd_errs, sd_err)
        end

        rates = log2.(sd_errs[1 : end-1] ./ sd_errs[2 : end])
        # println("Rates: ", rates)
        @test all(rates .> 1.9)
    end

    @testset "donating_region" begin
        U = u"m/s"
        L = u"m"
        T = L / U

        dt = .1T

        # Trivial case
        s = Segment((0, 0), (0, 1))
        u(x, y) = U * [1., 0.]
        drp, drm = donating_region(s, u, dt)
        @test drp == Quadrangle((0, 0), (0, 1), (-dt / T, 1), (-dt / T, 0))
        @test isnothing(drm)

        # Fix the reference volume (but no adjustment needed)
        α_ref = .1L^2
        drp, drm = donating_region(s, u, dt, α=α_ref)
        @test smeasure(drp) == α_ref
        @test drp == Pentagon((0, 0), (0, 1), (-dt / T, 1), (-dt / T, .5), (-dt / T, 0))
        @test isnothing(drm)

        # Fix the reference volume (adjustment needed)
        α_ref = .125L^2
        drp, drm = donating_region(s, u, dt, α=α_ref)
        @test smeasure(drp) == α_ref
        @test drp ≈ Pentagon((0, 0), (0, 1), (-dt / T, 1), (-1.5dt / T, .5), (-dt / T, 0))
        @test isnothing(drm)

        # Slightly less trivial case
        s = Segment((0, 0), (0, 1))
        u(x, y) = U * [1., 1.]
        drp, drm = donating_region(s, u, dt)
        @test drp == Quadrangle((0, 0), (0, 1), (-dt / T, 1 - dt / T), (-dt / T, -dt / T))
        @test isnothing(drm)

        # Fix the reference volume (but no adjustment needed)
        α_ref = .1L^2
        drp, drm = donating_region(s, u, dt, α=α_ref)
        @test smeasure(drp) ≈ α_ref
        @test drp == Pentagon((0, 0), (0, 1), (-dt / T, 1 - dt / T), ((-dt / T, .5 - dt / T)), (-dt / T, -dt / T))
        @test isnothing(drm)

        # Fix the reference volume (adjustment needed)
        α_ref = .125L^2
        drp, drm = donating_region(s, u, dt, α=α_ref)
        @test smeasure(drp) == α_ref
        @test isnothing(drm)

        # Nontrivial case
        u(x, y) = U * [-.1 + sin(x/L) * cos(3y/L), sin((x + y) / L)]
        s = Segment((0, 0), (.3, .2))
        drp, drm = donating_region(s, u, dt)
        @test isnothing(drp)
        @test isa(drm, Quadrangle)

        α_ref = 6E-3L^2
        drp, drm = donating_region(s, u, dt, α=α_ref)
        @test isapprox(smeasure(drp) + smeasure(drm), α_ref, rtol=10eps())

        # Trivial, but opposite sign
        s = Segment((0, 0), (0, 1))
        u(x, y) = U * [-1., 0.]

        # Fix the reference volume (but no adjustment needed)
        α_ref = -.1L^2
        drp, drm = donating_region(s, u, dt, α=α_ref)
        @test smeasure(drm) == α_ref
        @test drm == Pentagon((0, 0), (0, 1), (dt / T, 1), (dt / T, .5), (dt / T, 0))
        @test isnothing(drp)

        # Fix the reference volume (adjustment needed)
        α_ref = -.125L^2
        drp, drm = donating_region(s, u, dt, α=α_ref)
        @test smeasure(drm) == α_ref
        @test drm ≈ Pentagon((0, 0), (0, 1), (dt / T, 1), (1.5dt / T, .5), (dt / T, 0))
        @test isnothing(drp)

        # Self-intersection case
        s = Segment((0, 0), (0, 1))
        u(x, y) = U * [y/L-0.5, 0.]
        drp, drm = donating_region(s, u, dt)
        @test smeasure(drp) == 0.0125L^2
        @test smeasure(drm) == -0.0125L^2

        # Fix the reference volume (but no adjustment needed)
        α_ref = 0L^2
        drp, drm = donating_region(s, u, dt, α=α_ref)
        @test smeasure(drp) == 0.0125L^2
        @test smeasure(drm) == -0.0125L^2

        # Fix the reference volume (adjustment needed)
        α_ref = 0.01L^2
        drp, drm = donating_region(s, u, dt, α=α_ref)
        @test smeasure(drp) + smeasure(drm) ≈ α_ref
    end

end
