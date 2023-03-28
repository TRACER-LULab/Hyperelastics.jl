module Hyperelastics

using Reexport
@reexport using NonlinearContinua
@reexport using ADTypes
using InverseLangevinApproximations
# using LossFunctions
# using Optimization
using Tullio
using SpecialFunctions
# using DataInterpolations
using QuadGK
using ComponentArrays, LabelledArrays, StructArrays
using LinearAlgebra, Statistics
using Term

export HyperelasticUniaxialTest, HyperelasticBiaxialTest
export HyperelasticProblem
export predict
export parameters, parameter_bounds
export available_models

abstract type AbstractHyperelasticTest <: NonlinearContinua.AbstractMaterialTest end
abstract type AbstractHyperelasticModel <: NonlinearContinua.AbstractMaterialModel end
abstract type AbstractDataDrivenHyperelasticModel <: AbstractHyperelasticModel end
abstract type AbstractHyperelasticProblem end

struct InvariantForm end

function HyperelasticProblem end

include("invariants.jl")
include("material_tests.jl")
include("model_functions.jl")
# include("stress_functions.jl")
include("datasets.jl")

include("general_form_incompressible_isotropic_models.jl")
include("isotropic_incompressible_models.jl")
include("isotropic_compressible_models.jl")

include("data_driven.jl")
include("macro_micro_macro_model.jl")
include("average_chain_behavior.jl")



function available_models()
    exclude = [:HorganMurphy, :KhiemItskov, :GeneralCompressible, :LogarithmicCompressible, :GeneralMooneyRivlin]
    ns = filter(x -> x ∉ [:citation, :update_history, :update_history!], names(Hyperelastics))
    hyperelastic_models = filter(x -> typeof(getfield(Hyperelastics, x)) <: Union{DataType,UnionAll}, ns)
    hyperelastic_models = filter(x -> !(getfield(Hyperelastics, x) <: Hyperelastics.AbstractDataDrivenHyperelasticModel) && (getfield(Hyperelastics, x) <: Hyperelastics.AbstractHyperelasticModel), hyperelastic_models)
    hyperelastic_models_sym = filter(x -> !(x in exclude), hyperelastic_models)
    return hyperelastic_models_sym
end

end
