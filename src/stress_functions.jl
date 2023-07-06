
"""
`StrainEnergyDensity(ψ, λ⃗, p)`

Returns the strain energy density for the hyperelastic model `ψ` with the principle stretches `λ⃗` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `λ⃗`: Vector of principal stretches
- `p`: Model parameters
"""
function NonlinearContinua.StrainEnergyDensity(ψ::AbstractHyperelasticModel, λ⃗::Vector{T}, p) where {T}
    @error "$(typeof(ψ)) does not have a Strain Energy Density implemented"
end

"""
`StrainEnergyDensity(ψ, F, p)`

Returns a function for the strain energy density function for the hyperelastic model based on calculating the principal stretches of the deformation gradient, `F`. The eigen values are found by the following procedure:
```math
C = F^T \\cdot F
a = transpose(eigvecs(C))
C^\\ast = (U^\\ast)^2 = a^T \\cdot C \\cdot a
\\vec{\\lambda} = diag(U)
```

Fields:
- `ψ`: Hyperelastic model
- `F`: Deformation gradient matrix
- `p`: Model parameters
"""
function NonlinearContinua.StrainEnergyDensity(ψ::AbstractHyperelasticModel, F::Matrix{T}, p) where {T}
    C = transpose(F) * F
    λ⃗² = eigvals(C)
    λ⃗ = sqrt.(abs.(λ⃗²))
    return StrainEnergyDensity(ψ, λ⃗::Vector{T}, p)
end

"""
`StrainEnergyDensity(ψ, I⃗, p, InvariantForm())`

Returns the strain energy density for the model `ψ` with invariants `I⃗` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `I⃗`: vector of principal stretches or stretch invariants, respectively.
- `p`: parameters
- `InvariantForm()`
"""
function NonlinearContinua.StrainEnergyDensity(ψ::AbstractHyperelasticModel{T}, ::Vector{R}, p) where {T<:InvariantForm, R}
    @error "$(typeof(ψ)) does not have a stretch Invariant Form of Strain Energy Density implemented"
end

"""
`StrainEnergyDensity(ψ, F, p)`

Returns the strain energy density for the model `ψ` with deformation gradient `F` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `F`: Deformation gradient tensor
- `p`: parameters
- `InvariantForm()`
"""
function NonlinearContinua.StrainEnergyDensity(ψ::AbstractHyperelasticModel{T}, F::Matrix{R}, p) where {T <: InvariantForm, R}
    StrainEnergyDensity(ψ, [I₁(F), I₂(F), I₃(F)], p)
end

"""
`SecondPiolaKirchoffStressTensor(ψ::AbstractHyperelasticModel, λ⃗::AbstractVector, p; ad_type=AD.ForwardDiffBackend())`

Returns the second PK stress tensor for the hyperelastic model `ψ` with the principle stretches `λ⃗` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `λ⃗`: Vector of principal stretches
- `p`: Model parameters
"""
function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::AbstractHyperelasticModel, State, p; ad_type=nothing, kwargs...)
    _SecondPiolaKirchoffStressTensor(ψ, State, p, ad_type; kwargs...)
end

function _SecondPiolaKirchoffStressTensor(ψ, λ⃗, p, ad_type; kwargs...)
    @error "Please load a support derivative backend:
            1. ForwardDiff.jl
            2. Zygote.jl
            3. Enzyme.jl
            4. FiniteDiff.jl
            5. Or define a custom method for the material model"
end

"""
`SecondPiolaKirchoffStressTensor(ψ::AbstractHyperelasticModel, F::AbstractMatrix, p; adb=AD.ForwardDiffBackend())`

Returns the second PK stress tensor for the hyperelastic model `ψ` with the deformation gradient `F` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `F`: Deformation gradient tensor
- `p`: Model parameters
"""
function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, F::Matrix{T}, p; kwargs...) where {T}
    σ = CauchyStressTensor(ψ, F, p; kwargs...)
    S = sqrt(det(F' * F)) * inv(F) * σ
    return S
end

"""
`CauchyStressTensor(ψ::AbstractHyperelasticModel, λ⃗::AbstractVector, p; adb=AD.ForwardDiffBackend())`

Returns the Cauchy stress tensor for the hyperelastic model `ψ` with the principle stretches `λ⃗` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `λ⃗`: Vector of principal stretches
- `p`: Model parameters
- `adb`: Differentiation backend from `AbstractDifferentiation.jl`
    """
# NonlinearContinua.CauchyStressTensor(ψ::AbstractHyperelasticModel, State, p; ad_type=nothing, kwargs...) = _CauchyStressTensor(ψ, State, p, ad_type; kwargs...)

function NonlinearContinua.CauchyStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, λ⃗::Vector{T}, p; kwargs...) where {T}
    S = SecondPiolaKirchoffStressTensor(ψ, λ⃗, p; kwargs...)
    σ = S .* λ⃗ ./ J(λ⃗)
    return σ
end
# function _CauchyStressTensor(ψ::AbstractHyperelasticModel, λ⃗::Vector{T}, p, ad_type; kwargs...) where T
#     @error "Please load a support derivative backend:
#             1. ForwardDiff.jl
#             2. Zygote.jl
#             3. Enzyme.jl
#             4. FiniteDiff.jl
#             5. Or define a custom method for the material model"
# end

"""
`CauchyStressTensor(ψ::AbstractHyperelasticModel, F::AbstractMatrix, p; adb=AD.ForwardDiffBackend())`

Returns the Cauchy stress tensor for the hyperelastic model `ψ` with the deformation gradient `F` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `F`: Deformation gradient tensor
- `p`: Model parameters
- `adb`: Differentiation backend from `AbstractDifferentiation.jl`
"""
function NonlinearContinua.CauchyStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, F::Matrix{T}, p; kwargs...) where {T}
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
