module GeometricVOF

using FastGaussQuadrature: gausslegendre
using LinearAlgebra
using Meshes
using Optim
using Roots
using Unitful
using StaticArrays

export PlanarHS, LVIRA
export measure, intersect,  complement, shift, reconstruct, distance, donating_region,
    smeasure, symmetric_difference, normal, tangent

include("initialization.jl")
include("halfspace.jl")
include("reconstruction.jl")
include("donating_region.jl")
include("utils.jl")

end
