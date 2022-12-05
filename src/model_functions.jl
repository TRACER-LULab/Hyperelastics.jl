
"""
`StrainEnergyDensity(ψ, λ⃗, p)`

Returns the strain energy density for the hyperelastic model `ψ` with the principle stretches `λ⃗` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `λ⃗`: Vector of principal stretches
- `p`: Model parameters
"""
function NonlinearContinua.StrainEnergyDensity(ψ::AbstractHyperelasticModel, λ⃗::AbstractVector, p)
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
function NonlinearContinua.StrainEnergyDensity(ψ::AbstractHyperelasticModel, F::AbstractMatrix, p)
    C = transpose(F) * F
    λ⃗² = eigvals(C)
    λ⃗ = sqrt.(abs.(λ⃗²))
    return StrainEnergyDensity(ψ::AbstractHyperelasticModel, λ⃗, p)
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
function NonlinearContinua.StrainEnergyDensity(ψ::AbstractHyperelasticModel, I⃗::AbstractVector, p, I::InvariantForm)
    @error "$(typeof(ψ)) does not have a stretch Invariant Form of Strain Energy Density implemented"
end

"""
`StrainEnergyDensity(ψ, F, p, InvariantForm())`

Returns the strain energy density for the model `ψ` with deformation gradient `F` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `F`: Deformation gradient tensor
- `p`: parameters
- `InvariantForm()`
"""
function NonlinearContinua.StrainEnergyDensity(ψ::AbstractHyperelasticModel, F::AbstractMatrix, p, I::InvariantForm)
    NonlinearContinua.StrainEnergyDensity(ψ, [I₁(F), I₂(F), I₃(F)], p, I)
end

"""
`SecondPiolaKirchoffStressTensor(ψ::AbstractHyperelasticModel, λ⃗::AbstractVector, p; adb=AD.ForwardDiffBackend())`

Returns the second PK stress tensor for the hyperelastic model `ψ` with the principle stretches `λ⃗` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `λ⃗`: Vector of principal stretches
- `p`: Model parameters
- `adb`: Differentiation backend from `AbstractDifferentiation.jl`
"""
function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::AbstractHyperelasticModel, λ⃗::AbstractVector, p; adb=AD.ForwardDiffBackend())
    W(λ⃗) = StrainEnergyDensity(ψ, λ⃗, p)
    ∂W∂λ = AD.gradient(adb, W, λ⃗)[1]
    return ∂W∂λ
end

"""
`SecondPiolaKirchoffStressTensor(ψ::AbstractHyperelasticModel, F::AbstractMatrix, p; adb=AD.ForwardDiffBackend())`

Returns the second PK stress tensor for the hyperelastic model `ψ` with the deformation gradient `F` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `F`: Deformation gradient tensor
- `p`: Model parameters
- `adb`: Differentiation backend from `AbstractDifferentiation.jl`
"""
function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::AbstractHyperelasticModel, F::AbstractMatrix, p; adb=AD.ForwardDiffBackend())
    σ = CauchyStressTensor(ψ, F, p)
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
function NonlinearContinua.CauchyStressTensor(ψ::AbstractHyperelasticModel, λ⃗::AbstractVector, p; adb=AD.ForwardDiffBackend())
    W(λ⃗) = StrainEnergyDensity(ψ, λ⃗, p)
    ∂W∂λ = AD.gradient(adb, W, λ⃗)[1]
    σ = ∂W∂λ .* λ⃗ ./ J(λ⃗)
    return σ
end

"""
`CauchyStressTensor(ψ::AbstractHyperelasticModel, F::AbstractMatrix, p; adb=AD.ForwardDiffBackend())`

Returns the Cauchy stress tensor for the hyperelastic model `ψ` with the deformation gradient `F` with parameters `p`.

Fields:
- `ψ`: Hyperelastic model
- `F`: Deformation gradient tensor
- `p`: Model parameters
- `adb`: Differentiation backend from `AbstractDifferentiation.jl`
"""
function NonlinearContinua.CauchyStressTensor(ψ::AbstractHyperelasticModel, F::AbstractMatrix, p; adb=AD.ForwardDiffBackend())
    B = F * F'
    a = eigvecs(B)'
    B_prin = Diagonal(a * B * a')
    V_prin = sqrt.(B_prin)
    V = a' * V_prin * a
    R = inv(V) * F
    λ⃗ = sqrt.(diag(B_prin))
    W = StrainEnergyDensity(ψ, λ⃗, p)
    σ̂ = CauchyStressTensor(ψ, λ⃗, p) |> Diagonal
    σ = R * a' * σ̂ * a * R'
    return σ
end


"""
`parameters(ψ::AbstractHyperelasticModel)`

Returns a tuple of the parameters required for the model

Fields
- `ψ`: Hyperelastics model
"""
function parameters(ψ::AbstractHyperelasticModel)
    @error "$(typeof(ψ)) does not have a parameters function implemented"
end

"""
`parameter_bounds(ψ::AbstractHyperelasticModel, test::AbstractHyperelasticTest)`
`parameter_bounds(ψ::AbstractHyperelasticModel, tests::Vector{AbstractHyperelasticTest})`

Returns a tuple of the parameter bounds provided the experimental data and model

Fields
- `ψ`: Hyperelastic model
- `test` or `tests`: The test or vector of tests to use in finding the parameter bounds.
"""
function parameter_bounds(ψ::AbstractHyperelasticModel, test::AbstractHyperelasticTest)
    lb = nothing
    ub = nothing
    return (lb=lb, ub=ub)
end

function parameter_bounds(ψ::AbstractHyperelasticModel, tests::Vector{<:AbstractHyperelasticTest})
    bounds = map(Base.Fix1(parameter_bounds, ψ), tests)
    lbs = getfield.(bounds, :lb)
    ubs = getfield.(bounds, :ub)
    if !(eltype(lbs) <: Nothing)
        lb_ps = fieldnames(eltype(lbs))
        lb = map(p -> p .=> maximum(getfield.(lbs, p)), lb_ps) |> NamedTuple
    else
        lb = nothing
    end

    if !(eltype(ubs) <: Nothing)
        ub_ps = fieldnames(eltype(ubs))
        ub = map(p -> p .=> minimum(getfield.(ubs, p)), ub_ps) |> NamedTuple
    else
        ub = nothing
    end
    return (lb=lb, ub=ub)
end

"""
Returns a constraint equation for models were parameters bounds are interdependent
"""
function constraints(ψ::AbstractHyperelasticModel, data::AbstractHyperelasticTest)
    return nothing
end
