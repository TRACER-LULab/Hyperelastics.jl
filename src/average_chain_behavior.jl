export DataDrivenAverageChainBehavior

"""
Average Chain Behavior - Data-Driven

Adapted from the code provided in the article's supplementary material

Parameters:
`data`: Hyperelastic Data Object
`fchain(λch, pch)`: A function for interpolating the stress-stretch behavior of the chain

> Amores VJ, Benítez JM, Montáns FJ. Average-chain behavior of isotropic incompressible polymers obtained from macroscopic experimental data. A simple structure-based WYPiWYG model in Julia language. Advances in Engineering Software. 2019 Apr 1;130:41-57.
"""
struct DataDrivenAverageChainBehavior{T,S} <: AbstractDataDrivenHyperelasticModel
    data::T
    pchain::S
end

pchain(x, y) = BSplineInterpolation(y, x, 4, :Uniform, :Uniform)

function DataDrivenAverageChainBehavior(data::AbstractHyperelasticTest; pchain=pchain)
    λchain(λ⃗) = sqrt(I₁(λ⃗)/3)
    λ₁ = getindex.(data.data.λ, 1)
    λ₃ = getindex.(data.data.λ, 3)
    s₁ = getindex.(data.data.s, 1)
    pch = @. (s₁ * λ₁) / (λ₁^2 - λ₃^2) * λchain(data.data.λ⃗)
    pch[1] = pch[2]
    pchain = pchain(λchain.(collect.(data.data.λ⃗)), pch)
    DataDrivenAverageChainBehavior{typeof(data), typeof(pchain)}(data, pchain)
end

function NonlinearContinua.StrainEnergyDensity(ψ::DataDrivenAverageChainBehavior, λ⃗::Vector, p)
    λchain = sqrt(I₁(λ⃗) / 3)
    wchain = quadgk(ψ.pnchain(x), 1.0, λchain)[1]
    W = 3*wchain
    return W
end

function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::DataDrivenAverageChainBehavior, λ⃗::Vector, p)

end

function NonlinearContinua.CauchyStressTensor(ψ::DataDrivenAverageChainBehavior, λ⃗::Vector, p)

end
function
    parameters(::DataDrivenAverageChainBehavior)
    return nothing
end
