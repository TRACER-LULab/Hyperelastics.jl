"""
Generic Compressible Model

Parameters: ψ, κ

Model: ``\\psi(\\vec{\\lambda})+\\frac{\\kappa}{2}(J-1)^2``

Example Implementation: `StrainEnergyDensityFunction(Compressible(NeoHookean()), λ⃗, (κ=κ, μ=μ))`
"""

struct GeneralCompressible <: AbstractHyperelasticModel
    incompressible_model::AbstractHyperelasticModel
end

function StrainEnergyDensityFunction(ψ::GeneralCompressible, λ⃗, p)
    StrainEnergyDensityFunction(ψ.incompressible_model, λ⃗, p) + p.κ/2*(prod(λ⃗)-1)^2
end
