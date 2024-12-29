module GeometricVOF

using FastGaussQuadrature: gausslegendre
using LinearAlgebra
using Meshes: Point, Vec, Ngon, Triangle, Quadrangle, measure, to, vertices
using Roots
using Unitful
using StaticArrays

export Point, Vec, Ngon, Triangle, Quadrangle, measure, to, Plane

include("initialization.jl")
include("halfspace.jl")

end
