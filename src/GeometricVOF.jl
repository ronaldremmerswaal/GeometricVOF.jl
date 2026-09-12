module GeometricVOF

using FastGaussQuadrature: gausslegendre
using LinearAlgebra
using Meshes
import Meshes: normal
using Roots
using Optim
using Unitful
using StaticArrays
using Printf

export PlanarHS, Parabola, StaticNgon, StaticParabolicNgon, capacity
export measure, intersect, intersect!, complement, shift, reconstruct, reconstruct!, distance, donating_region,
    smeasure, moments, symmetric_difference, normal, normal!, tangent, tangent!, sorted_unique_approx, donating_region!,
    mof, pmof, plvira, prost

include("halfspace.jl")
include("initialization.jl")
include("parabola.jl")
include("reconstruction.jl")
include("donating_region.jl")
include("utils.jl")
include("optimization.jl")

end
