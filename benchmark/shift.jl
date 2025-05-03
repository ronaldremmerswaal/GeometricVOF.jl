using BenchmarkTools
using GeometricVOF
using Unitful
using StaticArrays
using Meshes

𝛈 = GeometricVOF.angle_to_normal(0.7)
quad = Quadrangle((-.1, -.1), (.9, -.1), (.9, .9), (-.1, .9))

workspace = StaticNgon(tri)
αvol = smeasure(quad) / π

shifts = MVector{30, Float64}(undef)

@benchmark shift($quad, $𝛈, $αvol; workspace=$workspace, shift_workspace=$shifts)
