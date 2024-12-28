using GeometricVOF
using Meshes
using Unitful
using Test

@testset "GeometricVOF.jl" begin
    @testset "levelset ∩ ngon" begin
        # Linear
        t = Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        Φ(x, y) = x - 0.5u"m"
        @test GeometricVOF.measure(Φ, t) == (.5 - 1/8) * u"m^2"
    end

end
