module NonlinearContinua

using LinearAlgebra
using Accessors

abstract type AbstractMaterialModel end
abstract type AbstractMaterialState end
abstract type AbstractMaterialTest end

export I₁, I₂, I₃, I1, I2, I3, J
export MaterialHistory

## Material Tests
"""
Predicts the model behavior for provided experimental test.
"""
function predict(ψ::AbstractMaterialModel, test::AbstractMaterialTest, ps)
    @error "Method not implemented for model $(typeof(ψ)) and test $(typeof(test))"
end
function predict(ψ::AbstractMaterialModel, tests::Vector{<:AbstractMaterialTest}, ps)
    f(test) = predict(ψ, test, ps)
    results = map(f, tests)
    return results
end
export predict
## Material Properties
export MaterialHistory, update_history, update_history!
struct MaterialHistory{T,S} <: AbstractMaterialState
    value::Vector{T}
    time::Vector{S}
    function MaterialHistory(values, times)
        new{eltype(values),eltype(times)}(values, times)
    end
end
value(history::MaterialHistory) = history.value
time(history::MaterialHistory) = history.time

function update_history!(history::MaterialHistory, value, time)
    push!(history.value, value)
    push!(history.time, time)
    return nothing
end

function update_history(history::MaterialHistory, value, time)
    history = @set history.value = vcat(history.value, [value])
    history = @set history.time = vcat(history.time, [time])
    return history
end

## Energy Models
for Model ∈ [
    :StrainEnergyDensity,
    :StrainEnergyDensity!,
]
    @eval export $Model
    @eval @inline function $Model(M::AbstractMaterialModel, S::AbstractMaterialState, P) end
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
    @eval @inline function $Tensor(M::AbstractMaterialModel, S::AbstractMaterialState, P) end
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
    @eval @inline function $Tensor(M::AbstractMaterialModel, S::AbstractMaterialState, P) end
end


## Strain Tensors
for Tensor ∈ [
    :GreenStrainTensor,
    :AlmansiStrainTensor,
    :GreenStrainTensor!,
    :AlmansiStrainTensor!,
]
    @eval export $Tensor
    @eval @inline function $Tensor(M::AbstractMaterialModel, S::AbstractMaterialState, P) end
end

## Time Dependent Tensors
# Deformation
for Tensor ∈ [
    :VelocityGradientTensor,
    :VelocityGradientTensor!,
]
    @eval export $Tensor
    @eval @inline function $Tensor(M::AbstractMaterialModel, S::AbstractMaterialState, P) end
end

## Electric Field Tensors

## Charge Displacement Tensors

## Tensor Invariant Calculations
I₁(T::AbstractMatrix) = tr(T)
I₂(T::AbstractMatrix) = 1 / 2 * (tr(T)^2 - tr(T^2))
I₃(T::AbstractMatrix) = det(T)
J(T::AbstractMatrix) = sqrt(det(T))
const I1 = I₁
const I2 = I₂
const I3 = I₃

## Precompile
# using SnoopPrecompile

end
