using GeometricVOF
using Meshes
using Unitful
using Test

@testset "GeometricVOF.jl" begin
    @testset "levelset ∩ ngon" begin
        # Linear
        t = Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        Φl(x, y) = x - 0.5u"m"
        @test GeometricVOF.measure(Φl, t) == .5(1 - 1/4)u"m^2"

        # Parabolic
        q = Quadrangle((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
        for p = 1 : 5
            Φp(x, y) = y - (.5 + .25(x/u"m")^p)u"m"
            @test GeometricVOF.measure(Φp, q) == (1/2 + 1/(4(p+1)))u"m^2"
        end
    end

end
