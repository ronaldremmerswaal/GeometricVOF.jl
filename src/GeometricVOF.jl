module GeometricVOF

using FastGaussQuadrature: gausslegendre
using LinearAlgebra
using Meshes
using Optim
using Roots
using Unitful
using StaticArrays
using Polynomials

export PlanarHS, LVIRA
export measure, intersect, intersect!,  complement, shift, reconstruct, distance, donating_region,
    smeasure, symmetric_difference, normal, tangent, sorted_unique_approx

include("initialization.jl")
include("halfspace.jl")
include("reconstruction.jl")
include("donating_region.jl")
include("utils.jl")
include("optimization.jl")

end
