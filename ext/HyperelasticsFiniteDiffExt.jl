module HyperelasticsFiniteDiffExt

import FiniteDiff: finite_difference_gradient
import NonlinearContinua
using Hyperelastics
using LinearAlgebra
using ADTypes

"""
`SecondPiolaKirchoffStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, λ⃗::AbstractVector, p; adb=AD.ForwardDiffBackend())`

Returns the second PK stress tensor for the hyperelastic model `ψ` with the principle stretches `λ⃗` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `λ⃗`: Vector of principal stretches
- `p`: Model parameters
- `adb`: Differentiation backend from `AbstractDifferentiation.jl`
"""
function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, λ⃗::AbstractVector, p, ad_type::AutoFiniteDiff; return_type=eltype(λ⃗), kwargs...)
    W(λ⃗) = StrainEnergyDensity(ψ, λ⃗, p)
    ∂W∂λ = finite_difference_gradient(W, λ⃗, ad_type.fdtype, return_type)
    return ∂W∂λ
end

"""
`SecondPiolaKirchoffStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, F::AbstractMatrix, p; adb=AD.ForwardDiffBackend())`

Returns the second PK stress tensor for the hyperelastic model `ψ` with the deformation gradient `F` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `F`: Deformation gradient tensor
- `p`: Model parameters
- `adb`: Differentiation backend from `AbstractDifferentiation.jl`
"""
function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, F::AbstractMatrix, p, ad_type::AutoFiniteDiff; return_type=eltype(λ⃗), kwargs...)
    σ = CauchyStressTensor(ψ, F, p, ad_type=ad_type, return_type=return_type, kwargs...)
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
function NonlinearContinua.CauchyStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, λ⃗::AbstractVector, p, ad_type::AutoFiniteDiff; return_type=eltype(λ⃗), kwargs...)
    W(λ⃗) = StrainEnergyDensity(ψ, λ⃗, p)
    ∂W∂λ = finite_difference_gradient(W, λ⃗, ad_type.fdtype, return_type, kwargs...)
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
function NonlinearContinua.CauchyStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, F::AbstractMatrix, p, ad_type::AutoFiniteDiff; return_type=eltype(λ⃗), kwargs...)
    B = F * F'
    a = eigvecs(B)'
    B_prin = Diagonal(a * B * a')
    V_prin = sqrt.(B_prin)
    V = a' * V_prin * a
    R = inv(V) * F
    λ⃗ = sqrt.(diag(B_prin))
    σ̂ = CauchyStressTensor(ψ, λ⃗, p, ad_type=ad_type, return_type=return_type, kwargs...) |> Diagonal
    σ = R * a' * σ̂ * a * R'
    return σ
end

end
