module GeometricVOF

using FastGaussQuadrature: gausslegendre
using LinearAlgebra
using Meshes
using Roots
using Unitful
using StaticArrays

export measure, intersect, PlanarHS, complement

include("initialization.jl")
include("halfspace.jl")

end
