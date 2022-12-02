# Data Driven Modelling

## Introduction
Sometimes, the constituitive models used in Hyperelastic modelling are not sufficient to capture the full behavior of engineered materials. To this end, a small set of data-driven models are included. 

## Setup
Loading the required packages and the [`Treloar1944Uniaxial`](#) dataset and create a test set for prediction:

{cell outputshow = false resultshow = false}
```julia
using Hyperelastics
using CairoMakie, MakiePublication
using DataInterpolations
treloar = Treloar1944Uniaxial()
test = HyperelasticUniaxialTest(
            collect(range(extrema(getindex.(treloar.data.λ,1))..., length = 100)),
            name="test"
        )
```

## Sussman-Bathe Model
The [`SussmanBathe`](#) model constructor takes three arguments:
1. A [`HyperelasticUniaxialTest`](#)
2. An `Int` for the number of terms in the series
3. An interpolation function constructor, `f(σ, λ)`, which returns an interpolation function for predicting the Cauchy stress from the stretch. 

For this example, `CubicSpline` from [`DataInterpolations.jl`](https://github.com/PumasAI/DataInterpolations.jl) will be used and `k` will be varied.

Predicting the behavior of the model:
{cell}
```julia
f = Figure()
ax = Axis(f[1,1], xlabel = "Stretch", ylabel = "Stress [MPa]")
itp = CubicSpline
for k in 1:10
    ψ = SussmanBathe(treloar, interpolant = itp, k = k)
    pred = predict(ψ, test, [])
    λ̂₁ = getindex.(pred.data.λ,1)
    ŝ₁ = getindex.(pred.data.s,1)
    lines!(ax, λ̂₁, ŝ₁, label = "k = $(k), itp = $(split(string(itp),".")[2])")
end
λ₁ = getindex.(treloar.data.λ,1)
s₁ = getindex.(treloar.data.s,1)
scatter!(ax, λ₁, s₁, label = "Treloar", color = :black)
axislegend(position = :lt, nbanks = 2)
f
```

## Average Chain Behavior Model
The [`DataDrivenAverageChainBehavior`](#)  asdf.