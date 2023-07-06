module NonlinearContinua

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
