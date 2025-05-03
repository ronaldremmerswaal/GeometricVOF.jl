module GeometricVOF

using FastGaussQuadrature: gausslegendre
using LinearAlgebra
using Meshes
using Optim
using Roots
using Unitful
using StaticArrays
using ForwardMethods

export PlanarHS, LVIRA, StaticNgon
export measure, intersect, intersect!,  complement, shift, reconstruct, distance, donating_region,
    smeasure, symmetric_difference, normal, tangent, sorted_unique_approx

include("halfspace.jl")
include("initialization.jl")
include("reconstruction.jl")
include("donating_region.jl")
include("utils.jl")
include("optimization.jl")

end
