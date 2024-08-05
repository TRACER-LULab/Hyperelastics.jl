export GeneralCompressible, LogarithmicCompressible

struct GeneralCompressible{T} <: AbstractCompressibleModel{T}
    incompressible::AbstractHyperelasticModel{T}
    """
    $(SIGNATURES)
    Generic Compressible Model

    # Model:
    ```math
    \\psi_{compressible} = \\psi_{incompressible}(\\vec{\\lambda}_{incompressible})+\\kappa(J-1)^2
    ```

    # Arguments:
    - ψ
        - Incompressible model

    # Parameters:
    - ψ
        - Incompressible model parameters (see model selected for parameter names)
    - κ
    """
    GeneralCompressible(ψ::AbstractIncompressibleModel{T}) where {T} = new{T}(ψ)
end

# Definitions
## Vector Forms
### Principal Value Form
function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::GeneralCompressible{T},
    λ⃗::Vector{S},
    p;
    kwargs...,
) where {T<:PrincipalValueForm,S}
    ψ_vol = (p.κ / 2) * (prod(λ⃗) - 1)^2
    ψ_dev = StrainEnergyDensity(ψ.incompressible, λ⃗ ./ cbrt(prod(λ⃗)), p.ψ; kwargs...)
    return ψ_vol + ψ_dev
end

function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    ψ::GeneralCompressible{T},
    λ⃗::Vector{S},
    p;
    kwargs...,
) where {T<:PrincipalValueForm,S}
    s_vol = p.κ * (prod(λ⃗) - 1) * prod(λ⃗) ./ λ⃗
    s_dev = SecondPiolaKirchoffStressTensor(
        ψ.incompressible,
        λ⃗ ./ cbrt(prod(λ⃗)),
        p.ψ;
        kwargs...,
    )
    return s_vol + s_dev
end

function ContinuumMechanicsBase.CauchyStressTensor(
    ψ::GeneralCompressible{T},
    λ⃗::Vector{S},
    p;
    kwargs...,
) where {T<:PrincipalValueForm,S}
    σ_vol = p.κ * (prod(λ⃗) - 1) * prod(λ⃗)
    σ_dev = CauchyStressTensor(ψ.incompressible, λ⃗ ./ cbrt(prod(λ⃗)), p.ψ; kwargs...)
    return σ_vol .+ σ_dev
end

## Invariant Form
function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::GeneralCompressible{T},
    I⃗::Vector{S},
    p;
    kwargs...,
) where {T<:InvariantForm,S}
    J = sqrt(I⃗[3])
    StrainEnergyDensity(
        ψ.incompressible,
        [I⃗[1] / (J^(2 / 3)), I⃗[2] * J^(2 / 3), I⃗[3] / (J^2)],
        p.ψ;
        kwargs...,
    ) + p.κ / 2 * (J - 1)^2
end

## Tensor Forms
### Principal Value Form
function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::GeneralCompressible{T},
    F::Matrix{S},
    p;
    kwargs...,
) where {T<:PrincipalValueForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)
    StrainEnergyDensity(ψ.incompressible, F̄, p.ψ; kwargs...) + p.κ / 2 * (J - 1)^2
end

function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    ψ::GeneralCompressible{T},
    F::Matrix{S},
    p;
    kwargs...,
) where {T<:PrincipalValueForm,S}
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

function ContinuumMechanicsBase.CauchyStressTensor(
    ψ::GeneralCompressible{T},
    F::Matrix{S},
    p;
    kwargs...,
) where {T<:PrincipalValueForm,S}
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
function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::GeneralCompressible{T},
    F::Matrix{S},
    p;
    kwargs...,
) where {T<:InvariantForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)
    StrainEnergyDensity(ψ.incompressible, [I₁(F̄), I₂(F̄), I₃(F̄)], p.ψ; kwargs...) +
    p.κ / 2 * (J - 1)^2
end

function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    ψ::GeneralCompressible{T},
    F::Matrix{S},
    p;
    kwargs...,
) where {T<:InvariantForm,S}
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

function ContinuumMechanicsBase.CauchyStressTensor(
    ψ::GeneralCompressible{T},
    F::Matrix{S},
    p;
    kwargs...,
) where {T<:InvariantForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)

    # Deviatoric Portion
    σ_dev = CauchyStressTensor(ψ.incompressible, F̄, p.ψ; kwargs...)

    # Hydrostatic Portion
    ∂ψ∂I₃ = p.κ / 2 * (1 - inv(J))
    σ_vol = 2 * J * ∂ψ∂I₃ * I

    # Total Stress
    return σ_dev + σ_vol
end

struct LogarithmicCompressible{T} <: AbstractCompressibleModel{T}
    incompressible::AbstractHyperelasticModel

    """
    $(SIGNATURES)

    Logarithmic Compressible Model

    # Model:
    ```math
    \\psi_{compressible} = \\psi_{incompressible}(\\vec{\\lambda}_{incompressible})+\\kappa(J\\log{J} - J)
    ```

    # Arguments:
    - ψ
        - Incompressible model

    # Parameters:
    - ψ
        - See Selected hyperelastic model for the required parameters.
    - κ

    """
    LogarithmicCompressible(ψ::AbstractIncompressibleModel{T}) where {T} = new{T}(ψ)
end

# Strain Energy Density Definitions
## Vector Definitions
### Principal Value Form
function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::LogarithmicCompressible{T},
    λ⃗::Vector{S},
    p,
) where {T<:PrincipalValueForm,S}
    StrainEnergyDensity(ψ.incompressible, λ⃗ ./ cbrt(J(λ⃗)), p.ψ) +
    p.κ * (J(λ⃗) * log(J(λ⃗)) - J(λ⃗))
end

function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    ψ::LogarithmicCompressible{T},
    λ⃗::Vector{S},
    p;
    kwargs...,
) where {T<:PrincipalValueForm,S}
    J = prod(λ⃗)
    s_vol = p.κ * log(J) * (J ./ λ⃗)
    s_dev = SecondPiolaKirchoffStressTensor(ψ.incompressible, λ⃗ ./ cbrt(J), p.ψ; kwargs...)
    return s_vol + s_dev
end

function ContinuumMechanicsBase.CauchyStressTensor(
    ψ::LogarithmicCompressible{T},
    λ⃗::Vector{S},
    p;
    kwargs...,
) where {T<:PrincipalValueForm,S}
    J = prod(λ⃗)
    σ_vol = p.κ * log(J) * J
    σ_dev = CauchyStressTensor(ψ.incompressible, λ⃗ ./ cbrt(J), p.ψ; kwargs...)
    return σ_vol .+ σ_dev
end

### Invariant Form
function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::LogarithmicCompressible{T},
    I⃗::Vector{S},
    p,
) where {T<:InvariantForm,S}
    J = sqrt(I⃗[3])
    StrainEnergyDensity(
        ψ.incompressible,
        [I⃗[1] / (J^(2 / 3)), I⃗[2] * J^(2 / 3), I⃗[3] / (J^2)],
        p.ψ,
    ) + p.κ * (J * log(J) - J)
end

## Tensor Definitions
### Principal Value Form
function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::LogarithmicCompressible{T},
    F::Matrix{S},
    p,
) where {T<:PrincipalValueForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)
    StrainEnergyDensity(ψ.incompressible, F̄, p.ψ) + p.κ * (J * log(J) - J)
end

function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    ψ::LogarithmicCompressible{T},
    F::Matrix{S},
    p;
    kwargs...,
) where {T<:PrincipalValueForm,S}
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

function ContinuumMechanicsBase.CauchyStressTensor(
    ψ::LogarithmicCompressible{T},
    F::Matrix{S},
    p;
    kwargs...,
) where {T<:PrincipalValueForm,S}
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
function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::LogarithmicCompressible{T},
    F::Matrix{S},
    p,
) where {T<:InvariantForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)
    StrainEnergyDensity(ψ.incompressible, F̄, p.ψ) + p.κ * (J * log(J) - J)
end

function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    ψ::LogarithmicCompressible{T},
    F::Matrix{S},
    p;
    kwargs...,
) where {T<:InvariantForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)
    s_dev = SecondPiolaKirchoffStressTensor(ψ.incompressible, F̄, p.ψ; kwargs...)
    ∂ψ∂I3 = p.κ * log(J^2) / (4 * J)
    s_vol = 2 * J^2 * ∂ψ∂I3 * inv(F)
    return s_dev + s_vol
end

function ContinuumMechanicsBase.CauchyStressTensor(
    ψ::LogarithmicCompressible{T},
    F::Matrix{S},
    p;
    kwargs...,
) where {T<:InvariantForm,S}
    J = det(F)
    F̄ = F ./ cbrt(J)
    σ_dev = CauchyStressTensor(ψ.incompressible, F̄, p.ψ; kwargs...)
    ∂ψ∂I3 = p.κ * log(J^2) / (4 * J)
    σ_vol = 2 * J * ∂ψ∂I3 * I
    return σ_dev + σ_vol
end


struct CHIPFoam <: AbstractCompressibleModel
    Wg::Function
    """
    $(SIGNATURES)

    CHIPFoam Model

    # Model:
    Refer to: Lewis M. A robust, compressible, hyperelastic constitutive model for the mechanical response of foamed rubber. Technische Mechanik-European Journal of Engineering Mechanics. 2016;36(1-2):88-101.

    # Arguments:
    - isothermal
        - Boolean to determine if the model is isothermal or adiabatic

    # Parameters:
    - Ĝ
    - K̂
    - Jb
    - pg
    - C10
    - φ₀,
    - K
    - p₀
    - γ

    """
    function CHIPFoam(isothermal::Bool=true)
        if isothermal
            Wg(Jg, p₀, φ₀, γ) = p₀ * φ₀ * (Jg - log(Jg) - 1)
        else
            Wg(Jg, p₀, φ₀, γ) = p₀ * φ₀ * (Jg - 1 / (γ - 1) * (γ - Jg^(1 - γ)))
        end
        new(Wg)
    end
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::CHIPFoam,
    I::Vector,
    (; Ĝ, K̂, Jb, pg, C10, φ₀, K, p₀, γ)
)
    I1 = I[1]
    J = sqrt(I[3])
    p̃ = pg + C10 * (φ₀^(1 / 3) * (4J - 4 + 5φ₀) / (J - 1 + φ₀) - (4J - 1) * (4J + 1) / (3J^(4 / 3)))
    Jm = exp(-p̃ / K)
    J̄ = J / Jm
    Jg = Jm * (J̄ - 1 + φ₀) / (φ₀)
    f = (2J - 1) / (cbrt(J)) + (2 - 2J + φ₀) * cbrt((φ₀) / (J - (1 - φ₀)))

    cbrt_φ₀ = cbrt(φ₀)
    cbrt_φ₀2 = cbrt_φ₀^(2)
    dJ̄dJ_1 = (-C10 - 4C10 * (cbrt_φ₀2) + K * (cbrt_φ₀2)) / (K * exp((5.0C10 - pg - 5.0C10 * (cbrt_φ₀)) / K) * (cbrt_φ₀2))

    W_LB = Ĝ / 2 * (I1 - 3) + K̂ * (
        (Jb - 1) * (J - (Jb + 1) / 2) + ((J - Jb) >= 0) * ((J - 1)^2 / 2 - (Jb - 1) * (J - (Jb + 1) / 2))
    )

    W_D = C10 * (Jm * (I1 * f - 3 * (1 - φ₀)) - J * (I1 - 3) * (1 - φ₀) * dJ̄dJ_1)

    W_M = (1 - φ₀) * K * (Jm * log(Jm) - Jm + 1)
    W_g = ψ.Wg(Jg, p₀, φ₀, γ)
    return W_LB + +W_D + W_M + W_g
end

function ContinuumMechanicsBase.CauchyStressTensor(
    ψ::CHIPFoam,
    F::Matrix,
    p;
    ad_type=nothing,
    kwargs...
)
    J = det(F)
    F̄ = F ./ cbrt(J)
    C̄ = F̄' * F̄
    Ī₁ = tr(C̄)
    Ī₂ = zero(Ī₁)
    I₃ = J^2
    ∂W∂Ī₁, _, ∂W∂I₃ = ∂ψ(ψ, [Ī₁, Ī₂, I₃], p, ad_type)
    ∂W∂Ī₁
    ∂W∂J = ∂W∂I₃ * 2 * J

    B̄ = F̄ * F̄'

    σ = 2 / J * ∂W∂Ī₁ * (B̄ - (LinearAlgebra.I * (Ī₁ / 3))) + ∂W∂J * LinearAlgebra.I
    return σ
end

function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    ψ::CHIPFoam,
    F::Matrix,
    p;
    ad_type=nothing,
    kwargs...
)
    σ = CauchyStressTensor(ψ, F, p)
    S = det(F) * inv(F) * σ' * inv(F)'
end
