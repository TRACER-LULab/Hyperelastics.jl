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
Another data-driven model is providedin [`DataDrivenAverageChainBehavior`](#) where the chain characteristics are determined from data and extrapolated out to the continuum response. For this example, the default interpolation of a \(4^{\text{th}}\) order B-spline is used. 

{cell}
```julia
f = Figure()
ax = Axis(f[1,1], xlabel = "Stretch", ylabel = "Stress [MPa]")
fchain(x, y) = BSplineInterpolation(y, x, 3, :Uniform, :Uniform)
ψ = DataDrivenAverageChainBehavior(treloar, fchain= fchain)
pred = predict(ψ, test, [])
λ̂₁ = getindex.(pred.data.λ,1)
ŝ₁ = getindex.(pred.data.s,1)
lines!(ax, λ̂₁, ŝ₁, label = "DataDrivenAverageChainBehavior")
λ₁ = getindex.(treloar.data.λ,1)
s₁ = getindex.(treloar.data.s,1)
scatter!(ax, λ₁, s₁, label = "Treloar", color = :black)
axislegend(position = :lt)
f
``` 

The chain behavior can be accessed with:

{cell}
```julia
f = Figure()
ax = Axis(f[1,1], xlabel = L"\lambda_{ch}", ylabel = L"f(\lambda_{ch})")
lines!(ax,1.0:0.01:4.4, ψ.fchain.(1.0:0.01:4.4))
f
```

## Macro-Micro-Macro Model

Lastly, the Macro-Micro-Macro model is provided. The model determines the average chain behavior based on an average over the microsphere for chain orientations. 

{cell}
```julia
λ₁ = [1.04,3.1]
tests = Kawabata1981.(λ₁)
n1 = 100
p₀ = range(0.0, 0.7, length=100)|>collect
λ_max = maximum(maximum.(map(x->maximum.(x.data.λ),tests)))
λs = collect(range(0.001, 6.0, length=n1))
PChain(u) = BSplineApprox(u, λs, 3, 12, :Uniform, :Uniform)
weights = ones(sum(map(test->length(test.data.s), tests)))

# Set the model Parameters
optimizer=LBFGS()

# Create the Model
ψ = MacroMicroMacro(tests, PChain, p₀, optimizer = optimizer);
preds = predict(ψ, Kawabata1981.([1.6, 1.9,2.5]), [])
## Compare the Experimental and Predicted stresses
f = Figure()
ax = Makie.Axis(f[1,1], xlabel = "Stretch [-]", ylabel = "Stress [MPa]")
for (pred,test) in zip(preds, Kawabata1981.([1.6, 1.9, 2.5]))
    lines!(
        ax,
        getindex.(pred.data.λ, 2),
        getindex.(pred.data.s, 1),
        label = "Predicted"
    )
    scatter!(
        ax,
        getindex.(test.data.λ, 2),
        getindex.(test.data.s, 1),
        label = "Experimental"
    )
end
axislegend(position = :rb)
f
```