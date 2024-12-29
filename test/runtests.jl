using GeometricVOF
using Meshes
using Unitful
using Test

@testset "GeometricVOF.jl" begin
    @testset "levelset ∩ ngon" begin
        # Linear
        t = Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        for Φl ∈ [(x, y) -> x - 0.5u"m", (x, y) -> x ≤ 0.5u"m"]
            @test GeometricVOF.measure(Φl, t) == .5(1 - 1/4)u"m^2"
        end

        # Monomial
        q = Quadrangle((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
        for p = 1 : 5
            for Φp ∈ [(x, y) -> y - (.5 + .25(x/u"m")^p)u"m", (x, y) -> y ≤ (.5 + .25(x/u"m")^p)u"m"]
                @test GeometricVOF.measure(Φp, q) == (1/2 + 1/(4(p+1)))u"m^2"
            end
        end

        # Circle
        ref(t) = refine(t, TriSubdivision())
        mesh = ref(ref(ref(q)))
        R = .25u"m"
        Φcircle(x, y) = (x - .5123u"m")^2 + (y - .45u"m")^2 - R^2
        M = 0.0u"m^2"
        for c ∈ mesh
            M += GeometricVOF.measure(Φcircle, c)
        end

        @test M ≈ π * R^2
    end

end
