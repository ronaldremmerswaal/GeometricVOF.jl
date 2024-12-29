using GeometricVOF
using Meshes
using Unitful
using Test

@testset verbose=true "GeometricVOF.jl" begin
    @testset "initialize" begin
        # Trivial: all inside
        t = Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        for Φl ∈ [(x, y) -> x - 2u"m", (x, y) -> x ≤ 2u"m"]
            @test measure(Φl, t) == measure(t)
        end

        # Trivial: all outside
        for Φl ∈ [(x, y) -> x + 2u"m", (x, y) -> x ≤ -2u"m"]
            @test measure(Φl, t) == 0u"m^2"
        end

        # Linear
        for Φl ∈ [(x, y) -> x - 0.5u"m", (x, y) -> x ≤ 0.5u"m"]
            @test measure(Φl, t) == .5(1 - 1/4)u"m^2"
        end

        # Monomial
        q = Quadrangle((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
        for p = 1 : 5
            for Φp ∈ [(x, y) -> y - (.5 + .25(x/u"m")^p)u"m", (x, y) -> y ≤ (.5 + .25(x/u"m")^p)u"m"]
                @test measure(Φp, q) == (1/2 + 1/(4(p+1)))u"m^2"
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

        @test M ≈ π * R^2

        # Edge case: colinear to edge, but all inside
        for Φl ∈ [(x, y) -> y - (1u"m" - x), (x, y) -> y ≤ 1u"m" - x]
            @test measure(Φl, t) == measure(t)
        end

        # Edge case: colinear to edge, all outside
        for Φl ∈ [(x, y) -> -y + (1u"m" - x), (x, y) -> y > 1u"m" - x]
            @test measure(Φl, t) < eps()u"m^2"
        end
    end

    @testset "HalfSpace" begin
        t = Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))

        # Trivial: all inside
        tp = t ∩ PlanarHS([1, 0], 2.)
        @test tp == t

        # Trivial: all outside
        tp = t ∩ PlanarHS([-1, 0], -2.)
        @test isnothing(tp)

        # Triangle to Quadrangle
        tp = t ∩ PlanarHS([1, 0], .5)
        @test tp == Quadrangle((.5, 0), (.5, .5), (0, 1), (0, 0))

        # Edge case: colinear to edge, but all inside
        tp = t ∩ PlanarHS([1, 1] / √2, √.5)
        @test tp == t

        # Edge case: colinear to edge, but all outside
        tp = t ∩ PlanarHS(-[1, 1] / √2, -√.5)
        @test isnothing(tp)
    end

end
