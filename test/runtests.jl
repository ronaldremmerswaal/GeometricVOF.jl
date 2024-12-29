using GeometricVOF
using Meshes
using Unitful
using Test

@testset verbose=true "GeometricVOF.jl" begin
    @testset "initialize" begin
        # Trivial: all inside
        t = Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        for Φ ∈ [(x, y) -> x - 2u"m", (x, y) -> x ≤ 2u"m"]
            @test measure(Φ, t) == measure(t)
        end

        # Trivial: all outside
        for Φ ∈ [(x, y) -> x + 2u"m", (x, y) -> x ≤ -2u"m"]
            @test measure(Φ, t) == 0u"m^2"
        end

        # Linear
        for Φ ∈ [(x, y) -> x - 0.5u"m", (x, y) -> x ≤ 0.5u"m"]
            @test measure(Φ, t) == .5(1 - 1/4)u"m^2"
        end

        # Non-convex case
        nc = Ngon((0, 0), (1, 0), (.5, .5), (1, 1), (0, 1))
        for Φ ∈ [(x, y) -> x - .75u"m", (x, y) -> x ≤ .75u"m"]
            @test measure(Φ, nc) == (3/4 - 1/16)u"m^2"
        end

        # Monomial
        q = Quadrangle((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
        for p = 1 : 5
            for Φ ∈ [(x, y) -> y - (.5 + .25(x/u"m")^p)u"m", (x, y) -> y ≤ (.5 + .25(x/u"m")^p)u"m"]
                @test measure(Φ, q) == (1/2 + 1/(4(p+1)))u"m^2"
            end
        end

        # Circle
        ref(t) = refine(t, TriSubdivision())
        mesh = ref(ref(ref(q)))
        R = .25u"m"
        Φcircle(x, y) = (x - .5123u"m")^2 + (y - .45u"m")^2 - R^2
        M = 0.0u"m^2"
        for c ∈ mesh
            M += measure(Φcircle, c)
        end

        @test isapprox(M, π * R^2, rtol=10eps())

        # Edge case: colinear to edge, but all inside
        for Φl ∈ [(x, y) -> y - (1u"m" - x), (x, y) -> y ≤ 1u"m" - x]
            @test measure(Φl, t) == measure(t)
        end

        # Edge case: colinear to edge, all outside
        for Φl ∈ [(x, y) -> -y + (1u"m" - x), (x, y) -> y > 1u"m" - x]
            @test measure(Φl, t) < eps()u"m^2"
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

            # Triangle to Quadrangle
            tp = t ∩ PlanarHS([1, 0], .5u"m")
            @test tp == Quadrangle((.5, 0), (.5, .5), (0, 1), (0, 0))

            # Non-convex case
            nc = Ngon((0, 0), (1, 0), (.5, .5), (1, 1), (0, 1))
            ncp = nc ∩ PlanarHS([1, 0], .75u"m")
            @test ncp == Ngon((.75, 0), (.75, .25), (.5, .5), (.75, .75), (.75, 1), (0, 1), (0, 0))

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

                    α = measure(c, p)
                    s_comp = shift(c, 𝛈, α)
                    if α == 0u"m^2"
                        @test s_comp ≤ shift_min
                    elseif α == measure(c)
                        @test s_comp ≥ shift_max
                    else
                        @test isapprox(s_comp, s_ref, rtol=100eps())
                    end
                end
            end
        end
    end

end
