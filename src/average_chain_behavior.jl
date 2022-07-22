export DataDrivenAverageChainBehavior

"""
Average Chain Behavior - Data-Driven [^1]

Adapted from the code provided in the article's supplementary material

Parameters:
`data`: Hyperelastic Data Object
`fchain(λch, pch)`: A function for interpolating the stress-stretch behavior of the chain

[^1] > Amores VJ, Benítez JM, Montáns FJ. Average-chain behavior of isotropic incompressible polymers obtained from macroscopic experimental data. A simple structure-based WYPiWYG model in Julia language. Advances in Engineering Software. 2019 Apr 1;130:41-57.
"""
struct DataDrivenAverageChainBehavior end
function StrainEnergyDensityFunction(ψ::DataDrivenAverageChainBehavior, (; data, pchain))
    # Calculate the values for fchain and λchain from the experimental data
    @assert typeof(data) ∈ [UniaxialHyperelasticData, BiaxialHyperelasticData]
    λchain(λ⃗) = sqrt(I₁(λ⃗) / 3)

    λ₁ = getindex.(data.λ⃗, 1)
    λ₃ = getindex.(data.λ⃗, 3)
    s₁ = getindex.(data.s⃗, 1)
    pch = @. (s₁ * λ₁) / (λ₁^2 - λ₃^2) * λchain(data.λ⃗)
    pch[1] = pch[2]

    _pnchain = pchain(λchain.(collect.(data.λ⃗)), pch)

    wchain(λch) = quadgk(x -> _pnchain(x), 1.0, λch)[1]
    W(λ⃗) = 3*wchain(λchain(λ⃗))
    return W
end

function parameters(ψ::DataDrivenAverageChainBehavior)
    return (:data, :pchain)
end
