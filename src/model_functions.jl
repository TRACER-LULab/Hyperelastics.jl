
"""
strain_energy_density(ψ, λ⃗, p)

Returns a function for the strain energy density function for the hyperelastic model `ψ` with parameters `p`.
> ψ = strain_energy_density(Gent(), (μ = 10, Jₘ = 100.0))
"""
function ContinuumModels.StrainEnergyDensity(ψ::AbstractHyperelasticModel, λ⃗::AbstractVector, p)
    @error "$(typeof(ψ)) does not have a Strain Energy Density Function implemented"
end

"""
strain_energy_density(ψ, F, p)

Returns a function for the strain energy density function for the hyperelastic model based on calculating the principal stretches of the deformation gradient, `F`. The eigen values are found by the following procedure:
``C = F^T \\cdot F``
``a = transpose(eigvecs(C))``
``C^\\ast = (U^\\ast)^2 = a^T \\cdot C \\cdot a``
``\\vec{\\lambda} = diag(U)``

"""
function ContinuumModels.StrainEnergyDensity(ψ::AbstractHyperelasticModel, F::AbstractMatrix, p)
    C = transpose(F) * F
    a = eigvecs(C)'
    C_prin = Diagonal(a * C * a')
    U_prin = sqrt.(C_prin)
    λ⃗ = diag(U_prin)
    return strain_energy_density(ψ::AbstractHyperelasticModel, λ⃗, p)
end

"""
strain_energy_density(ψ, λ⃗, p, InvariantForm())

Returns a function for the strain energy density function for the hyperelastic model `ψ` with parameters `p`.
> ψ = strain_energy_density(Gent(), (μ = 10, Jₘ = 100.0), InvariantForm())
"""
function ContinuumModels.StrainEnergyDensity(ψ::AbstractHyperelasticModel, I⃗::AbstractVector, p, I::InvariantForm)
    @error "$(typeof(ψ)) does not have a stretch Invariant Form of Strain Energy Density Function implemented"
end

"""
strain_energy_density(ψ, F, p, InvariantForm())

Returns a function for the strain energy density function for the hyperelastic model `ψ` with parameters `p` given a deformation gradient, `F` where the invariants are calculated.
> ψ = strain_energy_density(Gent(), (μ = 10, Jₘ = 100.0), InvariantForm())
"""
function ContinuumModels.StrainEnergyDensity(ψ::AbstractHyperelasticModel, F::AbstractMatrix, p, I::InvariantForm)
    ContinuumModels.StrainEnergyDensity(ψ, [I₁(F), I₂(F), I₃(F)], p, I)
end

"""
nominal_stress_function(ψ, λ⃗, p; adb = AD.ForwardDiffBackend())

Return a function for the nominal (2nd Piola Kirchoff) Stress Function  for the hyperelastic model `ψ` with parameters `p` for principal stretchs, `\\vec{\\lambda}}`. Defaults to using ForwardDiff for calculating the gradient of the strain energy density function. Implementing a new method for a model requires adding a new function with the type of the model.

> s = nominal_stress_function(Gent(), [2.0, 2.0, 1/4.0], (μ = 10, Jₘ = 100.0))
"""
function ContinuumModels.SecondPiolaKirchoffStressTensor(ψ::AbstractHyperelasticModel, λ⃗::AbstractVector, p; adb=AD.ForwardDiffBackend())
    W(λ⃗) = StrainEnergyDensity(ψ, λ⃗, p)
    ∂W∂λ = AD.gradient(adb, W, λ⃗)[1]
    return ∂W∂λ
end

"""
Return a function for the nominal (2nd Piola Kirchoff) Stress Function  for the hyperelastic model `ψ` with parameters `p` for deformation gradient tensor, `F`. Defaults to using ForwardDiff for calculating the gradient of the strain energy density function. Implementing a new method for a model requires adding a new function with the type of the model.

> s = nominal_stress_function(Gent(), [2 -2 0; 1 1 0; 0 0 1], (μ = 10, Jₘ = 100.0))
"""
function ContinuumModels.SecondPiolaKirchoffStressTensor(ψ::AbstractHyperelasticModel, F::AbstractMatrix, p; adb=AD.ForwardDiffBackend())
    σ = CauchyStressTensor(ψ, F, p)
    S = sqrt(det(F'*F)) * inv(F) * σ
    return S
end

"""
true_stress(ψ, p; adb = AD.ForwardDiffBackend())

Return a function for the true (Cauchy) Stress Function  for the hyperelastic model `ψ` with parameters `p`.
> σ = true_stress(Gent(), (μ = 10, Jₘ = 100.0))
"""
function ContinuumModels.CauchyStressTensor(ψ::AbstractHyperelasticModel, λ⃗::AbstractVector, p; adb=AD.ForwardDiffBackend())
    W(λ⃗) = StrainEnergyDensity(ψ, λ⃗, p)
    ∂W∂λ = AD.gradient(adb, W, λ⃗)[1]
    σ = ∂W∂λ .* λ⃗ ./ J(λ⃗)
    return σ
end

"""
true_stress(ψ, F, p; adb = AD.ForwardDiffBackend())

Return a function for the true (Cauchy) Stress Function  for the hyperelastic model `ψ` with parameters `p` for deformation gradient tensor, `F`. Defaults to using ForwardDiff for calculating the gradient of the strain energy density function. Implementing a new method for a model requires adding a new function with the type of the model.

> s = true_stress(Gent(), [2 -2 0; 1 1 0; 0 0 1], (μ = 10, Jₘ = 100.0))
"""
function ContinuumModels.CauchyStressTensor(ψ::AbstractHyperelasticModel, F::AbstractMatrix, p; adb=AD.ForwardDiffBackend())
    B = F * F'
    a = eigvecs(B)'
    B_prin = Diagonal(a * B * a')
    V_prin = sqrt.(B_prin)
    V = a' * V_prin * a
    R = inv(V) * F
    λ⃗ = sqrt.(diag(B_prin))
    σ̂ = CauchyStressTensor(ψ, λ⃗, p) |> Diagonal
    σ = R * a' * σ̂ * a * R'
    return σ
end


"""
Returns a tuple of the parameters required for the model
"""
function parameters(ψ::AbstractHyperelasticModel)
    @error "$(typeof(ψ)) does not have a parameters function implemented"
end

"""
Returns a tuple of the parameter bounds provided the experimental data and model
"""
function parameter_bounds(ψ::AbstractHyperelasticModel, data::AbstractHyperelasticData)
    lb = nothing
    ub = nothing
    return (lb=lb, ub=ub)
end

"""
Returns a constraint equation for models were parameters bounds are interdependent
"""
function constraints(ψ::AbstractHyperelasticModel, data::AbstractHyperelasticData)
    return nothing
end

"""
returns the bibtex entry for the model
"""
function citation(ψ::AbstractHyperelasticModel)
    @error "model $(ψ) does not have a citation defined in the bib file. Please file an issue to have the citation added"
end

function get_citation(s::AbstractString)
    bib_file = import_bibtex("../CITATIONS.bib")
    export_bibtex(bib_file[s])
end
