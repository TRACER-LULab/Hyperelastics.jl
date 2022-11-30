module Hyperelastics

using Reexport
@reexport using NonlinearContinua
using InverseLangevinApproximations
using LossFunctions, Optimization
using AbstractDifferentiation, ForwardDiff
using Tullio
using SpecialFunctions
using DataInterpolations
using QuadGK
using ComponentArrays, LabelledArrays, StructArrays
using LinearAlgebra, Statistics
using Term

export HyperelasticUniaxialTest, HyperelasticBiaxialTest, HyperelasticProblem
export predict
export citation, parameters, parameter_bounds

abstract type AbstractHyperelasticTest <: NonlinearContinua.AbstractMaterialTest end
abstract type AbstractHyperelasticModel <: NonlinearContinua.AbstractMaterialModel end
abstract type AbstractDataDrivenHyperelasticModel <: AbstractHyperelasticModel end
abstract type AbstractHyperelasticProblem end

struct InvariantForm end

include("data_types.jl")
include("invariants.jl")
include("material_tests.jl")
include("model_functions.jl")
include("datasets.jl")
include("general_form_incompressible_isotropic_models.jl")
include("isotropic_incompressible_models.jl")
include("isotropic_compressible_models.jl")
include("data_driven.jl")
include("macro_micro_macro_model.jl")
include("average_chain_behavior.jl")

end
