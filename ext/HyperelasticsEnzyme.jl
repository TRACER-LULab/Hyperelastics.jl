module HyperelasticsEnzyme

import NonlinearContinua
using Enzyme
using Hyperelastics
using LinearAlgebra

"""
`SecondPiolaKirchoffStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, λ⃗::AbstractVector, p; adb=AD.ForwardDiffBackend())`

Returns the second PK stress tensor for the hyperelastic model `ψ` with the principle stretches `λ⃗` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `λ⃗`: Vector of principal stretches
- `p`: Model parameters
- `mode`: Differentiation mode for `Enzyme`. Default = `Forward`
"""
function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, λ⃗::AbstractVector, p; mode=Forward, kwargs...)
    W(λ⃗) = StrainEnergyDensity(ψ, λ⃗, p)
    ∂W∂λ = gradient(mode, W, λ⃗, kwargs...)
    return ∂W∂λ
end

"""
`SecondPiolaKirchoffStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, F::AbstractMatrix, p; adb=AD.ForwardDiffBackend())`

Returns the second PK stress tensor for the hyperelastic model `ψ` with the deformation gradient `F` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `F`: Deformation gradient tensor
- `p`: Model parameters
- `mode`: Differentiation mode for `Enzyme`. Default = `Forward`
"""
function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, F::AbstractMatrix, p; mode = Forward, kwargs...)
    σ = CauchyStressTensor(ψ, F, p, mode=mode, kwargs...)
    S = sqrt(det(F' * F)) * inv(F) * σ
    return S
end

"""
`CauchyStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, λ⃗::AbstractVector, p; adb=AD.ForwardDiffBackend())`

Returns the Cauchy stress tensor for the hyperelastic model `ψ` with the principle stretches `λ⃗` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `λ⃗`: Vector of principal stretches
- `p`: Model parameters
- `adb`: Differentiation backend from `AbstractDifferentiation.jl`
    """
function NonlinearContinua.CauchyStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, λ⃗::AbstractVector, p; mode = Forward, kwargs...)
    W(λ⃗) = StrainEnergyDensity(ψ, λ⃗, p)
    ∂W∂λ = gradient(mode, W, λ⃗, kwargs...)
    σ = ∂W∂λ .* λ⃗ ./ J(λ⃗)
    return σ
end

"""
`CauchyStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, F::AbstractMatrix, p; adb=AD.ForwardDiffBackend())`

Returns the Cauchy stress tensor for the hyperelastic model `ψ` with the deformation gradient `F` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `F`: Deformation gradient tensor
- `p`: Model parameters
- `adb`: Differentiation backend from `AbstractDifferentiation.jl`
"""
function NonlinearContinua.CauchyStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, F::AbstractMatrix, p; mode = Forward, kwargs...)
    B = F * F'
    a = eigvecs(B)'
    B_prin = Diagonal(a * B * a')
    V_prin = sqrt.(B_prin)
    V = a' * V_prin * a
    R = inv(V) * F
    λ⃗ = sqrt.(diag(B_prin))
    W = StrainEnergyDensity(ψ, λ⃗, p)
    σ̂ = CauchyStressTensor(ψ, λ⃗, p, mode = mode, kwargs...) |> Diagonal
    σ = R * a' * σ̂ * a * R'
    return σ
end
end
