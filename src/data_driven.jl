export SussmanBathe
"""
`SussmanBathe(data::HyperelasticUniaxialTest;  interpolant=CubicSpline, k::Integer = 5)`

Model:
- See paper

Parameters:
- None

Fields:
- `data`: Hyperelastic Uniaxial test to be used for determining the interpolation
- `k`: Order of the summation in the model.
- `interpolant`: Function of the form, `f(s, λ)` which returns a function `f(λ) = s`

> Sussman T, Bathe KJ. A model of incompressible isotropic hyperelastic material behavior using spline interpolations of tension–compression test data. Communications in numerical methods in engineering. 2009 Jan;25(1):53-63.
"""
struct SussmanBathe{T,S} <: AbstractDataDrivenHyperelasticModel{PrincipalValueForm}
    w′::Function
    test::AbstractHyperelasticTest
    k::T
    itp::S
end

function SussmanBathe(data::HyperelasticUniaxialTest; interpolant=CubicSpline, k::Integer=5)
    σ̂ = interpolant(getindex.(data.data.s, 1) .* getindex.(data.data.λ, 1), getindex.(data.data.λ, 1))
    f(i, λ) = (σ̂(λ^((4.0)^(-i))) + σ̂(λ^(-0.5 * (4.0^(-i))))) / λ
    w′(λ) = sum(Base.Fix2(f, λ), 0:k)
    SussmanBathe{typeof(k),typeof(σ̂)}(w′, data, k, σ̂)
end

NonlinearContinua.StrainEnergyDensity(ψ::SussmanBathe, λ⃗::Vector{T}, p) where {T} = sum(x -> quadgk(ψ.w′, 1.0, x)[1], λ⃗)

function NonlinearContinua.StrainEnergyDensity(ψ::SussmanBathe, F::Matrix{T}, p)  where {T}
    λ⃗ = eigvals(F)
    return StrainEnergyDensity(ψ, λ⃗, p)
end

NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::SussmanBathe, λ⃗::Vector{T}, p; kwargs...) where {T} = ψ.w′.(λ⃗)

NonlinearContinua.CauchyStressTensor(ψ::SussmanBathe, λ⃗::Vector{T}, p; kwargs...) where {T} = ψ.w′.(λ⃗) .* λ⃗

parameters(ψ::SussmanBathe) = ()
