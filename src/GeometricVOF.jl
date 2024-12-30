module GeometricVOF

using FastGaussQuadrature: gausslegendre
using LinearAlgebra
using Meshes
using Optim
using Roots
using Unitful
using StaticArrays

export measure, intersect, PlanarHS, complement, shift, reconstruct, distance, LVIRA

include("initialization.jl")
include("halfspace.jl")
include("reconstruction.jl")

end
