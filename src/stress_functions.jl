"""
$(TYPEDSIGNATURES)

Returns the strain energy density for the hyperelastic model `ψ` with the principle stretches `λ⃗` with parameters `p`.

# Arguments:
- `ψ::AbstractHyperelasticModel`: Hyperelastic model
- `λ⃗::Vector`: Vector of principal stretches
- `p`: Model parameters
"""
function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::AbstractHyperelasticModel{T},
    ::Vector{R},
    p,
) where {T,R}
    return throw(
        ArgumentError("$(typeof(ψ)) does not have a Strain Energy Density implemented"),
    )
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::AbstractHyperelasticModel{T},
    ::Vector{R},
    p,
) where {T<:InvariantForm,R}
    return throw(
        ArgumentError(
            "$(typeof(ψ)) does not have a stretch Invariant Form of Strain Energy Density implemented",
        ),
    )
end

"""
$(TYPEDSIGNATURES)

Returns a function for the strain energy density function for the hyperelastic model based on calculating the principal stretches of the deformation gradient, `F`. The eigen values are found by the following procedure:
```math
C = F^T \\cdot F
a = transpose(eigvecs(C))
C^\\ast = (U^\\ast)^2 = a^T \\cdot C \\cdot a
\\vec{\\lambda} = diag(U)
```

# Arguments:
- `ψ::AbstractHyperelasticModel`: Hyperelastic model
- `F::Matrix`: Deformation gradient matrix
- `p`: Model parameters
"""
function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::AbstractHyperelasticModel{T},
    F::Matrix{R},
    p,
) where {T<:PrincipalValueForm,R}
    C = transpose(F) * F
    λ⃗² = eigvals(C)
    λ⃗ = sqrt.(abs.(λ⃗²))
    return StrainEnergyDensity(ψ, λ⃗::Vector{R}, p)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::AbstractHyperelasticModel{T},
    F::Matrix{R},
    p,
) where {T<:InvariantForm,R}
    StrainEnergyDensity(ψ, [I₁(F), I₂(F), I₃(F)], p)
end

"""
$(TYPEDSIGNATURES)

Returns the second PK stress tensor for the hyperelastic model `ψ` with the principle stretches `λ⃗` with parameters `p`.

# Arguments:
- `ψ::AbstractHyperelasticModel`: Hyperelastic model
- `λ⃗::Vector`: Vector of principal stretches
- `p`: Model parameters
- `ad_type`: Automatic differentiation backend (see `ADTypes.jl`)
"""
function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    ψ::AbstractHyperelasticModel{T},
    λ⃗::Vector{R},
    p;
    ad_type = nothing,
    kwargs...,
) where {T<:PrincipalValueForm,R}
    ∂ψ(ψ, λ⃗, p, ad_type; kwargs...)
end

"""
$(TYPEDSIGNATURES)

Returns the second PK stress tensor for the hyperelastic model `ψ` with the deformation gradient `F` with parameters `p`.

# Arguments:
- `ψ::AbstractHyperelasticModel`: Hyperelastic model
- `F::Matrix`: Deformation gradient matrix
- `p`: Model parameters
- `ad_type`: Automatic differentiation backend (see `ADTypes.jl`)
"""
function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    ψ::AbstractHyperelasticModel{T},
    F::Matrix{R},
    p;
    kwargs...,
) where {T<:PrincipalValueForm,R}
    σ = CauchyStressTensor(ψ, F, p; kwargs...)
    S = sqrt(det(F' * F)) * inv(F) * σ
    return S
end

function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    ψ::AbstractHyperelasticModel{T},
    F::Matrix{R},
    p;
    ad_type = nothing,
    kwargs...,
) where {T<:InvariantForm,R}
    C = F'*F
    I1 = I₁(C)
    I2 = I₂(C)
    I3 = I₃(C)
    ∂ψ∂I = ∂ψ(ψ, [I1, I2, I3], p, ad_type; kwargs...)
    S = 2∂ψ∂I[1] * F' + 2∂ψ∂I[2] * (I1 * F' + F' * F * F') + 2I3 * ∂ψ∂I[3] * inv(F)
    return S
end

function ∂ψ(ψ, λ⃗, p, ad_type; kwargs...)
    return throw(ArgumentError("Please load a support derivative backend:
            1. ForwardDiff.jl
            2. Zygote.jl
            3. Enzyme.jl
            4. FiniteDiff.jl
            5. Or define a custom method for the material model"))
end


"""
$(TYPEDSIGNATURES)
Returns the Cauchy stress tensor for the hyperelastic model `ψ` with the principle stretches `λ⃗` with parameters `p`.

# Arguments:
- `ψ::AbstractHyperelasticModel`: Hyperelastic model
- `λ⃗::Vector`: Vector of principal stretches
- `p`: Model parameters
- `ad_type`: Automatic differentiation backend (see `ADTypes.jl`)
"""
function ContinuumMechanicsBase.CauchyStressTensor(
    ψ::Hyperelastics.AbstractHyperelasticModel{T},
    λ⃗::Vector{R},
    p;
    kwargs...,
) where {T<:PrincipalValueForm,R}
    S = SecondPiolaKirchoffStressTensor(ψ, λ⃗::Vector{R}, p; kwargs...)
    σ = S .* λ⃗
    return σ
end

"""
$(TYPEDSIGNATURES)

Returns the Cauchy stress tensor for the hyperelastic model `ψ` with the deformation gradient `F` with parameters `p`.

# Arguments:
- `ψ::AbstractHyperelasticModel`: Hyperelastic model
- `F::Matrix`: Deformation gradient matrix
- `p`: Model parameters
- `ad_type`: Automatic differentiation backend (see `ADTypes.jl`)
"""
function ContinuumMechanicsBase.CauchyStressTensor(
    ψ::Hyperelastics.AbstractHyperelasticModel{T},
    F::Matrix{S},
    p;
    kwargs...,
) where {T<:PrincipalValueForm,S}
    B = F * F'
    a = eigvecs(B)'
    B_prin = Diagonal(a * B * a')
    V_prin = sqrt.(B_prin)
    V = a' * V_prin * a
    R = inv(V) * F
    λ⃗ = sqrt.(diag(B_prin))
    σ̂ = CauchyStressTensor(ψ, λ⃗, p; kwargs...) |> Diagonal
    σ = R * a' * σ̂ * a * R'
    return σ
end

function ContinuumMechanicsBase.CauchyStressTensor(
    ψ::Hyperelastics.AbstractHyperelasticModel{T},
    F::Matrix{S},
    p;
    ad_type,
    kwargs...,
) where {T<:InvariantForm,S}
    C = F'*F
    I1 = I₁(C)
    I2 = I₂(C)
    I3 = I₃(C)
    J = sqrt(I3)
    B = F * F'
    ∂ψ∂I = ∂ψ(ψ, [I1, I2, I3], p, ad_type; kwargs...)
    σ =
        2 * inv(J) * (∂ψ∂I[1] + I1 * ∂ψ∂I[2]) * B - 2 * inv(J) * ∂ψ∂I[2] * B^2 +
        2 * J * ∂ψ∂I[3] * I
    return σ
end

function ContinuumMechanicsBase.MaterialElasticStiffnessTensor(
 ψ::Hyperelastics.AbstractHyperelasticModel{T},
    F::Matrix{S},
    p;
    ad_type,
    kwargs...,
) where {T<:InvariantForm,S}
    I1 = I₁(F)
    I2 = I₂(F)
    I3 = I₃(F)
    J = sqrt(I3)
    B = F * F'
    ∂²ψ∂I² = ∂²ψ(ψ, [I1, I2, I3], p, ad_type; kwargs...)
    C = zeros(3, 3, 3, 3)
    for i in 1:3
        for j in 1:3
            for k in 1:3
                for l in 1:3
                    C[i, j, k, l] =
                        2 * inv(J) * ∂²ψ∂I²[1][i, j, k, l] +
                        2 * I1 * ∂²ψ∂I²[2][i, j, k, l] +
                        2 * I2 * ∂²ψ∂I²[3][i, j, k, l] +
                        2 * J * ∂²ψ∂I²[4][i, j, k, l]
                end
            end
        end
    end
    return C
end
