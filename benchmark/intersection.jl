using BenchmarkTools
using GeometricVOF
using Unitful
using StaticArrays
using Meshes

hs = PlanarHS{2}(SVector(1., 0.), 0.333u"m")
tri = Triangle((0., 0.), (1., 0.), (0., 1.))
quad = Quadrangle((-.1, -.1), (.9, -.1), (.9, .9), (-.1, .9))

memory = MVector{30, eltype(tri.vertices)}(undef)

@benchmark intersect($quad, $hs)
@benchmark intersect!($memory, $quad, $hs)
