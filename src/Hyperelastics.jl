module Hyperelastics

using InverseLangevinApproximations
using LossFunctions
using Optimization
using AbstractDifferentiation, ForwardDiff
using Tullio
using SpecialFunctions
using ComponentArrays
using SymbolicUtils
using Bibliography

export AbstractHyperelasticData, UniaxialHyperelasticData, BiaxialHyperelasticData, HyperelasticProblem
export I₁, I₂, I₃, J
export StrainEnergyDensityFunction, NominalStressFunction, TrueStressFunction, citation, parameters, parameter_bounds

abstract type AbstractHyperelasticData end
abstract type AbstractHyperelasticModel end

struct InvariantForm end
include("data_types.jl")
include("invariants.jl")
include("model_functions.jl")
include("datasets.jl")
include("isotropic_incompressible_models.jl")
include("data_driven.jl")
include("macro_micro_macro_model.jl")
include("average_chain_behavior.jl")
include("optimization_interface.jl")

end
