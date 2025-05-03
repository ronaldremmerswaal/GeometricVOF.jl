using BenchmarkTools
using GeometricVOF
using Unitful
using StaticArrays
using Meshes

hs = PlanarHS(SVector(1., 0.), 0.333u"m")
tri = Triangle((0., 0.), (1., 0.), (0., 1.))
quad = Quadrangle((-.1, -.1), (.9, -.1), (.9, .9), (-.1, .9))

out = StaticNgon(tri)
workspace = StaticNgon(tri)

@benchmark intersect($quad, $tri)
@benchmark intersect!($out, $quad, $tri; workspace=$workspace)
