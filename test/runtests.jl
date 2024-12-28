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
        Φp(x, y) = y - (.5 + .25(x/u"m")^2)u"m"
        @test GeometricVOF.measure(Φp, q) == (1/2 + 1/12)u"m^2"
    end

end
