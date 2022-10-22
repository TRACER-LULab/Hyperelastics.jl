module Hyperelastics

using InverseLangevinApproximations, ContinuumModels
using LossFunctions, Optimization
using AbstractDifferentiation, ForwardDiff
using Tullio
using SpecialFunctions
using DataInterpolations
using QuadGK
using ComponentArrays, LabelledArrays
using LinearAlgebra

export UniaxialHyperelasticData, BiaxialHyperelasticData, HyperelasticProblem
export citation, parameters, parameter_bounds

abstract type AbstractHyperelasticData end
abstract type AbstractHyperelasticModel <: ContinuumModels.AbstractMaterialModel end
abstract type AbstractDataDrivenHyperelasticModel <: AbstractHyperelasticModel end
abstract type AbstractHyperelasticProblem end

struct InvariantForm end

include("data_types.jl")
include("invariants.jl")
include("model_functions.jl")
include("datasets.jl")
include("isotropic_incompressible_models.jl")
include("isotropic_compressible_models.jl")
include("data_driven.jl")
include("macro_micro_macro_model.jl")
include("average_chain_behavior.jl")
include("optimization_interface.jl")

# include("precompile.jl")

end
