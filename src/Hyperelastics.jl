module Hyperelastics

using Reexport
@reexport using NonlinearContinua
@reexport using ADTypes
using InverseLangevinApproximations
using LossFunctions

# using Tullio
using SpecialFunctions
using DataInterpolations
using QuadGK
using ComponentArrays, LabelledArrays, StructArrays
using LinearAlgebra, Statistics
using Term

export HyperelasticUniaxialTest, HyperelasticBiaxialTest
export HyperelasticProblem
export predict
export parameters, parameter_bounds
export InvariantForm, PrincipalValueForm, DataDrivenForm

abstract type AbstractHyperelasticTest{T,S} <: NonlinearContinua.AbstractMaterialTest end

abstract type AbstractHyperelasticModel{T} <: NonlinearContinua.AbstractMaterialModel end

abstract type AbstractIncompressibleModel{T} <: AbstractHyperelasticModel{T} end
abstract type AbstractCompressibleModel{T} <: AbstractHyperelasticModel{T} end
abstract type AbstractDataDrivenHyperelasticModel{T} <: AbstractHyperelasticModel{T} end

abstract type AbstractHyperelasticProblem end

struct InvariantForm end
struct PrincipalValueForm end
struct DataDrivenForm end

"""
`HyperelasticProblem(ψ::AbstractHyperelasticModel, test::AbstractHyperelasticTest, u₀, ps=Nothing;
    adb=AD.ForwardDiffBackend(), loss=L2DistLoss(), adtype=Optimization.AutoForwardDiff())`

`HyperelasticProblem(ψ::AbstractHyperelasticModel, tests::Vector{<:AbstractHyperelasticTest}, u₀, ps=Nothing;
    adb=AD.ForwardDiffBackend(), loss=L2DistLoss(), adtype=Optimization.AutoForwardDiff())`

Creates an `OptimizationProblem` for use in [`Optimization.jl`](https://docs.sciml.ai/Optimization/stable/) to find the optimal parameters.

Fields:
- `ψ`: material model to use
- `test` or `tests`: A single or vector of hyperelastics tests to use when fitting the parameters
- `u₀`: Initial guess for parameters
- `ps`: Any additional parameters for calling predict
- `adb`: Select differentiation type from [`ADTypes.jl`](https://github.com/SciML/ADTypes.jl). The type is automatically applied to the type of AD applied to the Optimization Problem also.
- `loss`: Loss function from [`LossFunctions.jl`](https://github.com/JuliaML/LossFunctions.jl)
"""
struct HyperelasticProblem{iip,F,uType,P,LB,UB,I,LC,UC,S,K}
    f::F
    u0::uType
    p::P
    lb::LB
    ub::UB
    int::I
    lcons::LC
    ucons::UC
    sense::S
    kwargs::K
end

# include("../ext/HyperelasticsOptimizationExt.jl")
include("invariants.jl")
include("material_tests.jl")
include("model_functions.jl")
include("stress_functions.jl")
include("datasets.jl")

include("general_form_incompressible_isotropic_models.jl")
include("isotropic_incompressible_models.jl")
include("isotropic_compressible_models.jl")

include("data_driven.jl")
include("average_chain_behavior.jl")
# include("macro_micro_macro_model.jl")
end
