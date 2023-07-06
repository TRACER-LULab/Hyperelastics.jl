
module Hyperelastics

using Reexport
@reexport module NonlinearContinua

using LinearAlgebra
using RecursiveArrayTools

abstract type AbstractMaterialModel end
abstract type AbstractMaterialState end
abstract type AbstractMaterialTest end

export I₁, I₂, I₃, J
export MaterialHistory, update_history, update_history!
export predict

## Material Tests
"""
`predict(ψ::AbstractMaterialModel, test::AbstractMaterialTest, ps)`

Fields:
- `ψ`: Material Model
- `test` or `tests`: A single test or vector of tests. This is used to predict the response of the model in comparison to the experimental data provided.
- `ps`: Model parameters to be used.
"""
function predict(ψ::AbstractMaterialModel, test::AbstractMaterialTest, ps)
    @error "Method not implemented for model $(typeof(ψ)) and test $(typeof(test))"
end
function predict(ψ::AbstractMaterialModel, tests::Vector{<:AbstractMaterialTest}, ps, args...)
    f(test) = predict(ψ, test, ps,args...)
    results = map(f, tests)
    return results
end

"""
`MaterialHistory(values::Vector, times::Vector)`

Structure for storing the behavior of a material as it evolves in time. Design to be used in time-dependent models such as viscoelasticity.

"""
struct MaterialHistory{T} <: AbstractMaterialState
    value::VectorOfArray
    time::Vector{T}
    function MaterialHistory(value::Vector, time::T) where { T}
        new{T}(VectorOfArray([value]), [time])
    end
    function MaterialHistory(value::Matrix, time::T) where {T}
        new{T}(VectorOfArray([value]), [time])
    end
end

## Energy Models
for Model ∈ [
    :StrainEnergyDensity,
    :StrainEnergyDensity!,
]
    @eval export $Model
    @eval @inline function $Model(M::AbstractMaterialModel, S, P; kwargs...) end
end
## Stress Tensors
for Tensor ∈ [
    :FirstPiolaKirchoffStressTensor,
    :SecondPiolaKirchoffStressTensor,
    :CauchyStressTensor,
    :FirstPiolaKirchoffStressTensor!,
    :SecondPiolaKirchoffStressTensor!,
    :CauchyStressTensor!,
]
    @eval export $Tensor
    @eval @inline function $Tensor(M::AbstractMaterialModel, S, P, ;kwargs...) end
end

## Deformation Tensors
for Tensor ∈ [
    :DeformationGradientTensor,
    :InverseDeformationGradientTensor,
    :RightCauchyGreenDeformationTensor,
    :LeftCauchyGreenDeformationTensor,
    :InverseLeftCauchyGreenDeformationTensor,
    :DeformationGradientTensor!,
    :InverseDeformationGradientTensor!,
    :RightCauchyGreenDeformationTensor!,
    :LeftCauchyGreenDeformationTensor!,
    :InverseLeftCauchyGreenDeformationTensor!,
]
    @eval export $Tensor
    @eval @inline function $Tensor(M::AbstractMaterialModel, S::AbstractMaterialState, P;kwargs...) end
end


## Strain Tensors
for Tensor ∈ [
    :GreenStrainTensor,
    :AlmansiStrainTensor,
    :GreenStrainTensor!,
    :AlmansiStrainTensor!,
]
    @eval export $Tensor
    @eval @inline function $Tensor(M::AbstractMaterialModel, S::AbstractMaterialState, P;kwargs...) end
end

## Time Dependent Tensors
# Deformation
for Tensor ∈ [
    :VelocityGradientTensor,
    :VelocityGradientTensor!,
]
    @eval export $Tensor
    @eval @inline function $Tensor(M::AbstractMaterialModel, S::AbstractMaterialState, P; kwargs...) end
end

## Electric Field Tensors

## Charge Displacement Tensors

## Tensor Invariant Calculations
I₁(T::AbstractMatrix) = tr(T)
I₂(T::AbstractMatrix) = 1 / 2 * (tr(T)^2 - tr(T^2))
I₃(T::AbstractMatrix) = det(T)
J(T::AbstractMatrix) = sqrt(det(T))

end

@reexport using ADTypes
@reexport module InverseLangevinApproximations

export CohenRounded3_2, CohenExact3_2, PusoApproximation, TreloarApproximation, WarnerApproximation, KuhnGrunApproximation, BergstromApproximation, PadeApproximation_1_4, PadeApproximation_1_2, PadeApproximation_5_0, PadeApproximation_3_0, Jedynak2017, ArrudaApproximation

CohenRounded3_2(y) = y * (3 - y^2) / (1 - y^2)

ArrudaApproximation(y) = 3y + 9 / 5 * y + 297 / 175 * y^3 + 297 / 175 * y^5 + 1539 / 875 * y^7 + 126117 / 67375 * y^9

PadeApproximation_3_2(y) = y * (3 - 36 / 35 * y^2) / (1 - 33 / 35 * y^2)

PusoApproximation(y) = 3 * y / (1 - y^3)

TreloarApproximation(y) = 3 * y / (1 - (3 / 5 * y^2 + 36 / 175 * y^4 + 108 / 875 * y^6))

WarnerApproximation(y) = 3 * y / (1 - y^2)

KuhnGrunApproximation(y) = 3y + 9y^3 / 5 + 297y^5 / 175 + 1539y^7 / 875 + 126117y^9 / 67375 + 43733439y^11 / 21896875 + 231321177y^13 / 109484375 + 20495009043y^15 / 9306171875 + 1073585186448381y^17 / 476522530859375 + 4387445039583y^19 / 1944989921875

function BergstromApproximation(y)
    if abs(y) < 0.84136
        return 1.31446tan(1.58986y) + 0.91209y
    elseif 0.84136 <= abs(y) < 1.0
        return 1 / (sign(y) - y)
    end
end

PadeApproximation_1_4(y) = 3y / (1 - 3y^2 / 5 - 36y^4 / 175)

PadeApproximation_1_2(y) = 3y / (1 - 3y^2 / 5)

PadeApproximation_5_0(y) = 3y + 9y^3 / 5 + 297y^5 / 175

PadeApproximation_3_0(y) = 3y + 9y^3 / 5

Jedynak2017(y) = y * (3 - 773 / 768 * y^2 - 1300 / 1351 * y^4 + 501 / 340 * y^6 - 678 / 1385 * y^8) / (1 - y) / (1 + 866 / 853 * y)

end
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
