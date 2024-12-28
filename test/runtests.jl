using GeometricVOF
using Meshes
using Unitful
using Test

@testset "GeometricVOF.jl" begin
    @testset "initialization" begin
        t = Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        Φ(x, y) = x - 0.5u"m"
        GeometricVOF.volume(Φ, t)
    end

end
