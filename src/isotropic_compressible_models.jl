export GeneralCompressible, LogarithmicCompressible
"""
Generic Compressible Model

Model:
```math
\\psi_{compressible} = \\psi_{incompressible}(\\vec{\\lambda}_{incompressible})+\\kappa(J-1)^2
```

Parameters:
- ψ
    - Incompressible model parameters (see model selected for parameter names)
- κ

"""
struct GeneralCompressible{T} <: AbstractCompressibleModel{T}
    incompressible::AbstractHyperelasticModel{T}
    GeneralCompressible(W::AbstractIncompressibleModel{T}) where {T} = new{T}(W)
end

# Definitions
## Vector Forms
### Principal Value Form
function NonlinearContinua.StrainEnergyDensity(ψ::GeneralCompressible{T}, λ⃗::Vector{S}, p; kwargs...) where {T<:PrincipalValueForm, S}
    ψ_vol = (p.κ / 2) * (prod(λ⃗) - 1)^2
    ψ_dev = StrainEnergyDensity(ψ.incompressible, λ⃗ ./ cbrt(prod(λ⃗)), p.ψ; kwargs...)
    return ψ_vol + ψ_dev
end

function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::GeneralCompressible{T}, λ⃗::Vector{S}, p; kwargs...) where {T<:PrincipalValueForm,S}
    s_vol = p.κ * (prod(λ⃗) - 1) * prod(λ⃗) ./ λ⃗
    s_dev = SecondPiolaKirchoffStressTensor(ψ.incompressible, λ⃗ ./ cbrt(prod(λ⃗)), p.ψ; kwargs...)
    return s_vol + s_dev
end

function NonlinearContinua.CauchyStressTensor(ψ::GeneralCompressible{T}, λ⃗::Vector{S}, p; kwargs...) where {T<:PrincipalValueForm,S}
    σ_vol = p.κ * (prod(λ⃗) - 1) * prod(λ⃗)
    σ_dev = CauchyStressTensor(ψ.incompressible, λ⃗ ./ cbrt(prod(λ⃗)), p.ψ; kwargs...)
    return σ_vol .+ σ_dev
end

## Invariant Form
function NonlinearContinua.StrainEnergyDensity(ψ::GeneralCompressible{T}, I⃗::Vector{S}, p; kwargs...) where {T<:InvariantForm,S}
    J = sqrt(I⃗[3])
    StrainEnergyDensity(ψ.incompressible, [I⃗[1]/(J^(2/3)), I⃗[2]*J^(2/3), I⃗[3]/(J^2)], p.ψ; kwargs...) + p.κ / 2 * (J - 1)^2
end

## Tensor Forms
### Principal Value Form
function NonlinearContinua.StrainEnergyDensity(ψ::GeneralCompressible{T}, F::Matrix{S}, p; kwargs...) where {T<:PrincipalValueForm, S}
    J = det(F)
    F̄ = F ./ cbrt(J)
    StrainEnergyDensity(ψ.incompressible, F̄, p.ψ; kwargs...) + p.κ / 2 * (J - 1)^2
end

function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::GeneralCompressible{T}, F::Matrix{S}, p; kwargs...) where {T<:PrincipalValueForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)

    # Deviatoric Portion
    s_dev = SecondPiolaKirchoffStressTensor(ψ.incompressible, F̄, p.ψ; kwargs...)

    # Hydrostatic Portion
    B = F * F'
    a = eigvecs(B)'
    B_prin = Diagonal(a * B * a')
    V_prin = sqrt.(B_prin)
    V = a' * V_prin * a
    R = inv(V) * F
    σ̂ = p.κ * (J - 1) * J
    σ = R * a' * σ̂ * I * a * R'
    s_vol = sqrt(det(F' * F)) * inv(F) * σ

    # Total Stress
    return s_dev .+ s_vol
end

function NonlinearContinua.CauchyStressTensor(ψ::GeneralCompressible{T}, F::Matrix{S}, p; kwargs...) where {T<:PrincipalValueForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)

    # Deviatoric Portion
    σ_dev = CauchyStressTensor(ψ.incompressible, F̄, p.ψ; kwargs...)

    # Hydrostatic Portion
    B = F * F'
    a = eigvecs(B)'
    B_prin = Diagonal(a * B * a')
    V_prin = sqrt.(B_prin)
    V = a' * V_prin * a
    R = inv(V) * F
    σ̂ = p.κ * (J - 1) * J
    σ_vol = R * a' * σ̂ * I * a * R'

    # Total Stress
    return σ_dev .+ σ_vol
end


### Invariant Form
function NonlinearContinua.StrainEnergyDensity(ψ::GeneralCompressible{T}, F::Matrix{S}, p; kwargs...) where {T<:InvariantForm, S}
    J = det(F)
    F̄ = F ./ cbrt(J)
    StrainEnergyDensity(ψ.incompressible, [I₁(F̄), I₂(F̄), I₃(F̄)], p.ψ; kwargs...) + p.κ / 2 * (J - 1)^2
end

function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::GeneralCompressible{T}, F::Matrix{S}, p; kwargs...) where {T<:InvariantForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)

    # Deviatoric Portion
    s_dev = SecondPiolaKirchoffStressTensor(ψ.incompressible, F̄, p.ψ; kwargs...)

    # Hydrostatic Portion
    ∂ψ∂I₃ = p.κ / 2 * (1 - inv(J))
    s_vol = 2 * J^2 * ∂ψ∂I₃ * inv(F)

    # Total Stress
    return s_dev .+ s_vol
end

function NonlinearContinua.CauchyStressTensor(ψ::GeneralCompressible{T}, F::Matrix{S}, p; kwargs...) where {T<:InvariantForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)

    # Deviatoric Portion
    σ_dev = CauchyStressTensor(ψ.incompressible, F̄, p.ψ; kwargs...)

    # Hydrostatic Portion
    ∂ψ∂I₃ = p.κ / 2 * (1 - inv(J))
    σ_vol = 2*J*∂ψ∂I₃*I

    # Total Stress
    return σ_dev + σ_vol
end

"""
Logarithmic Compressible Model

Model:
```math
\\psi_{compressible} = \\psi_{incompressible}(\\vec{\\lambda}_{incompressible})+\\kappa(J\\log{J} - J)
```

Parameters:
- ψ
    - See Selected hyperelastic model for the required parameters.
- κ

"""
struct LogarithmicCompressible{T} <: AbstractCompressibleModel{T}
    incompressible::AbstractHyperelasticModel
    LogarithmicCompressible(W::AbstractIncompressibleModel{T}) where {T} = new{T}(W)
end

# Strain Energy Density Definitions
## Vector Definitions
### Principal Value Form
function NonlinearContinua.StrainEnergyDensity(ψ::LogarithmicCompressible{T}, λ⃗::Vector{S}, p) where {T<:PrincipalValueForm, S}
    StrainEnergyDensity(ψ.incompressible, λ⃗./cbrt(J(λ⃗)), p.ψ) + p.κ * (J(λ⃗) * log(J(λ⃗)) - J(λ⃗))
end

function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::LogarithmicCompressible{T}, λ⃗::Vector{S}, p; kwargs...) where {T<:PrincipalValueForm,S}
    J = prod(λ⃗)
    s_vol = p.κ * log(J) * (J ./ λ⃗)
    s_dev = SecondPiolaKirchoffStressTensor(ψ.incompressible, λ⃗ ./ cbrt(J), p.ψ; kwargs...)
    return s_vol + s_dev
end

function NonlinearContinua.CauchyStressTensor(ψ::LogarithmicCompressible{T}, λ⃗::Vector{S}, p; kwargs...) where {T<:PrincipalValueForm,S}
    J = prod(λ⃗)
    σ_vol = p.κ * log(J) * J
    σ_dev = CauchyStressTensor(ψ.incompressible, λ⃗ ./ cbrt(J), p.ψ; kwargs...)
    return σ_vol .+ σ_dev
end

### Invariant Form
function NonlinearContinua.StrainEnergyDensity(ψ::LogarithmicCompressible{T}, I⃗::Vector{S}, p) where {T<:InvariantForm,S}
    J = sqrt(I⃗[3])
    StrainEnergyDensity(ψ.incompressible, [I⃗[1] / (J^(2 / 3)), I⃗[2] * J^(2 / 3), I⃗[3] / (J^2)], p.ψ) + p.κ * (J * log(J) - J)
end

## Tensor Definitions
### Principal Value Form
function NonlinearContinua.StrainEnergyDensity(ψ::LogarithmicCompressible{T}, F::Matrix{S}, p) where {T<:PrincipalValueForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)
    StrainEnergyDensity(ψ.incompressible, F̄, p.ψ) + p.κ * (J * log(J) - J)
end

function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::LogarithmicCompressible{T}, F::Matrix{S}, p; kwargs...) where {T<:PrincipalValueForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)

    # Deviatoric Portion
    s_dev = SecondPiolaKirchoffStressTensor(ψ.incompressible, F̄, p.ψ; kwargs...)

    # Hydrostatic Portion
    B = F * F'
    a = eigvecs(B)'
    B_prin = Diagonal(a * B * a')
    V_prin = sqrt.(B_prin)
    V = a' * V_prin * a
    R = inv(V) * F
    σ̂ = p.κ * log(J) * J
    σ = R * a' * σ̂ * I * a * R'
    s_vol = sqrt(det(F' * F)) * inv(F) * σ

    # Total Stress
    return s_dev .+ s_vol
end

function NonlinearContinua.CauchyStressTensor(ψ::LogarithmicCompressible{T}, F::Matrix{S}, p; kwargs...) where {T<:PrincipalValueForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)

    # Deviatoric Portion
    σ_dev = CauchyStressTensor(ψ.incompressible, F̄, p.ψ; kwargs...)

    # Hydrostatic Portion
    B = F * F'
    a = eigvecs(B)'
    B_prin = Diagonal(a * B * a')
    V_prin = sqrt.(B_prin)
    V = a' * V_prin * a
    R = inv(V) * F
    λ⃗ = sqrt.(diag(B_prin))
    σ̂ = p.κ * log(J) * J
    σ_vol = R * a' * σ̂ * I * a * R'

    # Total Stress
    return σ_dev .+ σ_vol
end

### Invariant Form
function NonlinearContinua.StrainEnergyDensity(ψ::LogarithmicCompressible{T}, F::Matrix{S}, p) where {T<:InvariantForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)
    StrainEnergyDensity(ψ.incompressible, F̄, p.ψ) + p.κ * (J * log(J) - J)
end

function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::LogarithmicCompressible{T}, F::Matrix{S}, p; kwargs...) where {T<:InvariantForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)
    s_dev = SecondPiolaKirchoffStressTensor(ψ.incompressible, F̄, p.ψ; kwargs...)
    ∂ψ∂I3 = p.κ*log(J^2)/(4*J)
    s_vol = 2*J^2*∂ψ∂I3*inv(F)
    return s_dev + s_vol
end

function NonlinearContinua.CauchyStressTensor(ψ::LogarithmicCompressible{T}, F::Matrix{S}, p; kwargs...) where {T<:InvariantForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)
    σ_dev = CauchyStressTensor(ψ.incompressible, F̄, p.ψ; kwargs...)
    ∂ψ∂I3 = p.κ*log(J^2)/(4*J)
    σ_vol = 2*J*∂ψ∂I3*I
    return σ_dev + σ_vol
end
