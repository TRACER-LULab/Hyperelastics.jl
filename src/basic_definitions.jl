
struct InvariantForm end

"""
First stretch invariant - Currently requires the addition of 5 times the machine precision to allow AD to work correctly

``I_1(\\vec{\\lambda}) = \\lambda_1^2+\\lambda_2^2+\\lambda_3^2 + 5\\varepsilon``
"""
function I₁(λ⃗)
    sum(λ⃗ .^ 2) + 10eps(Float64)
end

"""
Second Stretch invariant

``I_2(\\vec{\\lambda}) = \\lambda_1^{-2}+\\lambda_2^{-2}+\\lambda_3^{-2}``
"""
function I₂(λ⃗)
    sum(λ⃗ .^ (-2))
end

"""
Third Stretch invariant

``I_3(\\vec{\\lambda}) = (\\lambda_1\\lambda_2\\lambda_3)^2``
"""
function I₃(λ⃗)
    prod(λ⃗)^2
end

"""
Volumetric Stretch

``J(\\vec{\\lambda}) = \\lambda_1\\lambda_2\\lambda_3``
"""
function J(λ⃗)
    prod(λ⃗)
end

"""
StrainEnergyDensityFunction(ψ, p)

Returns a function for the strain energy density function for the hyperelastic model `ψ` with parameters `p`.
> W = StrainEnergyDensityFunction(Gent(), (μ = 10, Jₘ = 100.0))
"""
function StrainEnergyDensityFunction(ψ, p)
    @error "$(typeof(ψ)) does not have a Strain Energy Density Function implemented"
end

"""
StrainEnergyDensityFunction(ψ, p)

Returns a function for the strain energy density function for the hyperelastic model `ψ` with parameters `p`.
> W = StrainEnergyDensityFunction(Gent(), (μ = 10, Jₘ = 100.0))
"""
function StrainEnergyDensityFunction(ψ, p, I::InvariantForm)
    @error "$(typeof(ψ)) does not have a stretch Invariant Form of Strain Energy Density Function implemented"
end

"""
NominalStressFunction(ψ, p; adb = AD.ForwardDiffBackend())

Return a function for the nominal (2nd Piola Kirchoff) Stress Function  for the hyperelastic model `ψ` with parameters `p`. Defaults to using ForwardDiff for calculating the gradient of the strain energy density function. Implementing a new method for a model requires adding a new function with the type of the model.
> s = NominalStressFunction(Gent(), (μ = 10, Jₘ = 100.0))
"""
function NominalStressFunction(ψ, p; adb=AD.ForwardDiffBackend())
    W = StrainEnergyDensityFunction(ψ, p)
    function s(λ⃗)
        ∂W∂λᵢ = AD.gradient(adb, W, λ⃗)[1]
        s⃗ᵢ = ∂W∂λᵢ .- ∂W∂λᵢ[3] .* λ⃗[3] ./ λ⃗
        return s⃗ᵢ
    end
    return s
end

"""
TrueStressFunction(ψ, p; adb = AD.ForwardDiffBackend())

Return a function for the true (Cauchy) Stress Function  for the hyperelastic model `ψ` with parameters `p`.
> σ = TrueStressFunction(Gent(), (μ = 10, Jₘ = 100.0))
"""
function TrueStressFunction(ψ, p; adb=AD.ForwardDiffBackend())
    W = StrainEnergyDensityFunction(ψ, p)
    function σ(λ⃗)
        ∂W∂λᵢ = map(λ⃗ᵢ -> AD.gradient(adb, W, λ⃗ᵢ)[1], λ⃗)
        s⃗ᵢ = ∂W∂λᵢ .- getindex.(∂W∂λᵢ, 3) .* getindex.(λ⃗, 3) .* map(λ⃗ᵢ -> 1 ./ λ⃗ᵢ, λ⃗)
        σ⃗ᵢ = map(x -> x[1] .* x[2], zip(s⃗ᵢ, λ⃗))
        return σ⃗ᵢ
    end
    return σ
end

function parameters(ψ)
    @error "$(typeof(ψ)) does not have a parameters function implemented"
end
