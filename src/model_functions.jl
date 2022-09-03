
"""
StrainEnergyDensityFunction(ψ, p)

Returns a function for the strain energy density function for the hyperelastic model `ψ` with parameters `p`.
> ψ = StrainEnergyDensityFunction(Gent(), (μ = 10, Jₘ = 100.0))
"""
function StrainEnergyDensityFunction(ψ::AbstractHyperelasticModel, λ⃗, p)
    @error "$(typeof(ψ)) does not have a Strain Energy Density Function implemented"
end

"""
StrainEnergyDensityFunction(ψ, p, InvariantForm())

Returns a function for the strain energy density function for the hyperelastic model `ψ` with parameters `p`.
> ψ = StrainEnergyDensityFunction(Gent(), (μ = 10, Jₘ = 100.0), InvariantForm())
"""
function StrainEnergyDensityFunction(ψ::AbstractHyperelasticModel, I⃗, p, I::InvariantForm)
    @error "$(typeof(ψ)) does not have a stretch Invariant Form of Strain Energy Density Function implemented"
end

"""
NominalStressFunction(ψ, p; adb = AD.ForwardDiffBackend())

Return a function for the nominal (2nd Piola Kirchoff) Stress Function  for the hyperelastic model `ψ` with parameters `p`. Defaults to using ForwardDiff for calculating the gradient of the strain energy density function. Implementing a new method for a model requires adding a new function with the type of the model.

> s = NominalStressFunction(Gent(), (μ = 10, Jₘ = 100.0))
"""
function NominalStressFunction(ψ::AbstractHyperelasticModel, λ⃗, p; adb=AD.ForwardDiffBackend())
    W(λ⃗) = StrainEnergyDensityFunction(ψ, λ⃗, p)
    ∂W∂λ = AD.gradient(adb, W, λ⃗)[1]
    return ∂W∂λ
    # @tullio Δs[i,j] := (∂W∂λ[i]*λ⃗[i] - ∂W∂λ[j]*λ⃗[j])/λ⃗[i]
    # # sᵢ = ∂W∂λ .- ∂W∂λ[3] .* λ⃗[3] ./ λ⃗
    # # return sᵢ
    # Δs
end

"""
TrueStressFunction(ψ, p; adb = AD.ForwardDiffBackend())

Return a function for the true (Cauchy) Stress Function  for the hyperelastic model `ψ` with parameters `p`.
> σ = TrueStressFunction(Gent(), (μ = 10, Jₘ = 100.0))
"""
function TrueStressFunction(ψ::AbstractHyperelasticModel, λ⃗, p; adb=AD.ForwardDiffBackend())
    W(λ⃗) = StrainEnergyDensityFunction(ψ, λ⃗, p)
    ∂W∂λ = AD.gradient(adb, W, λ⃗)[1]
    σ = ∂W∂λ.*λ⃗
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
    return (lb = lb, ub = ub)
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
