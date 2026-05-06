module SemiDiscretizationMethod

using Reexport
@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using StaticArrays
@reexport using Arpack

using ForwardDiff
using KrylovKit
using QuadGK
using Lazy: iterated, take

include("structures_method.jl")
include("structures_input.jl")
include("structures_result.jl")

include("functions_utility.jl")
include("functions_discretization.jl")
include("functions_method.jl")

include("functions_LRmapping.jl")
include("sdm_periodic_sol.jl")

export SemiDiscretization, NumericSD, 
ProportionalMX,
Delay,DelayMX,
Additive,
LDDEProblem,
DiscreteMapping, DiscreteMapping_1step,
DiscreteMapping_LR,
fixPointOfMapping, spectralRadiusOfMapping,
extract_SDM_system,get_periodic_solution
end # module
