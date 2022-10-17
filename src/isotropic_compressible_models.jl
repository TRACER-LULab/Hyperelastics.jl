# export GeneralCompressible
# """
# Generic Compressible Model

# Parameters: ψ, κ

# Model: ``\\psi(\\vec{\\lambda})+\\kappa(J-1)^2``

# Example Implementation: `StrainEnergyDensityFunction(Compressible(NeoHookean()), λ⃗, (κ=κ, μ=μ))`
# """
# struct GeneralCompressible <: AbstractHyperelasticModel
#     incompressible_model::AbstractHyperelasticModel
# end

# function StrainEnergyDensityFunction(ψ::GeneralCompressible, λ⃗, p)
#     StrainEnergyDensityFunction(ψ.incompressible_model, λ⃗, p) + p.κ/2 * (prod(λ⃗) - 1)^2
# end

# function TrueStressFunction(ψ::GeneralCompressible, λ⃗, p; adb=AD.ForwardDiffBackend())
#     σ_dev =  p.κ *(prod(λ⃗) - 1)
#     σ_1 = TrueStressFunction(ψ.incompressible_model, λ⃗./cbrt(prod(λ⃗)), p).*λ⃗./cbrt(prod(λ⃗))./prod(λ⃗)
#     return σ_1 .+ σ_dev
# end


# σ = 152.9633929, 148.51830355, 148.51830355
# T1 = 2.9633929, -1.48169645, -1.48169645
