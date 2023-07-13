export DataDrivenAverageChainBehavior

"""
`DataDrivenAverageChainBehavior(data::HyperelasticUniaxialTest; fchain=fchain)`

Model:
- Adapted from the code provided in the article's supplementary material

Parameters:
- None
Fields:
- `data`: A Uniaxial Hyperelastic Test
- `fchain(λch, pch)`: A constructor for an approximation with the form f(x, y) => f̂(x) = y

> Amores VJ, Benítez JM, Montáns FJ. Average-chain behavior of isotropic incompressible polymers obtained from macroscopic experimental data. A simple structure-based WYPiWYG model in Julia language. Advances in Engineering Software. 2019 Apr 1;130:41-57.
> Amores VJ, Benítez JM, Montáns FJ. Data-driven, structure-based hyperelastic manifolds: A macro-micro-macro approach to reverse-engineer the chain behavior and perform efficient simulations of polymers. Computers & Structures. 2020 Apr 15;231:106209.
"""
struct DataDrivenAverageChainBehavior{T,S} <:
       AbstractDataDrivenHyperelasticModel{PrincipalValueForm}
    data::T
    fchain::S
end

fchain(x, y) = BSplineInterpolation(y, x, 3, :Uniform, :Uniform)

function DataDrivenAverageChainBehavior(data::HyperelasticUniaxialTest; fchain = fchain)
    λ₁ = getindex.(data.data.λ, 1)
    s₁ = getindex.(data.data.s, 1)
    f̂ = @. s₁ / (λ₁ - 1 / (λ₁^2))
    f̂[1] = f̂[2]
    λch = @. sqrt(I₁(data.data.λ) / 3)
    fch = fchain(λch, f̂)
    DataDrivenAverageChainBehavior{typeof(data.name),typeof(fch)}(data.name, fch)
end

function NonlinearContinua.StrainEnergyDensity(
    ψ::DataDrivenAverageChainBehavior,
    λ⃗::Vector,
    p,
)
    λchain = sqrt(I₁(λ⃗) / 3)
    W = ψ.fchain(λchain) / 2 * I₁(λ⃗)
    return W
end

function NonlinearContinua.SecondPiolaKirchoffStressTensor(
    ψ::DataDrivenAverageChainBehavior,
    λ⃗::Vector,
    p;
    kwargs...,
)
    λch = sqrt(I₁(λ⃗) / 3)
    sᵢ = λ⃗ .* ψ.fchain(λch)
    return sᵢ
end

function NonlinearContinua.CauchyStressTensor(
    ψ::DataDrivenAverageChainBehavior,
    λ⃗::Vector,
    p;
    kwargs...,
)
    λch = sqrt(I₁(λ⃗) / 3)
    σᵢ = λ⃗ .^ 2 .* ψ.fchain(λch)
    return σᵢ
end

function parameters(::DataDrivenAverageChainBehavior)
    return nothing
end
