# Bayesian Parameter Estimation of Model Parameters
{.subtitle}
Learn how to use Bayesian parameter estimation with Hyperelastics.jl

While the traditional loss-based types of analysis generally work for predicting the parameters of hyperelastic materials, other alternatives for parameter existmation exist. One such example is Bayesian parameter fitting which provides a probablity of model parameters and noise in the measurements. [Turing.jl](https://github.com/TuringLang/Turing.jl) is used in this example:

{cell}
```julia
using Hyperelastics
using Turing, Distributions, LinearAlgebra
using MCMCChains, SciMLExpectations, KernelDensity
using CairoMakie, MakiePublication
using ComponentArrays:ComponentVector
using Optimization, OptimizationOptimJL

Turing.setadbackend(:forwarddiff)
```

For this example, we will use all of [`Kawabata1981`](#) Biaxial data and the [`ExpLn`](#) model:
{cell}
```julia
kawabata_data = map(λ₁ -> Kawabata1981(λ₁), [1.040, 1.060, 1.080, 1.100, 1.120, 1.14, 1.16, 1.2, 1.24, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4, 3.7])
ψ = ExpLn()
```

First a guess of the optimal parameters is performed with the `Optimization.jl` to find a starting guess for the mean values of the distributions:

{cell}
```julia
p₀ = ComponentVector(A = 0.198881173072011, a = 0.009823455836135064, b = 0.09660741942475941)
prob = HyperelasticProblem(ψ, kawabata_data, p₀)
sol = solve(prob, LBFGS())
```

{cell}
```julia
sol.retcode
```
 
To ease the setting of model parameters, a parametric struct is created to be passed to the Hyperelastic model:

{cell}
```julia
struct p_ExpLn{S, T, R}
    A::S
    a::T
    b::R
end
```

A helper function is created to get the stress predictions for a single data instance:

{cell}
```julia
function ŝ(ψ, test, p)
    pred = Hyperelastics.predict(ψ, test, p)
    s = getindex.(pred.data.s, 1)
    return s
end
```

Next the model is created:

{cell}
```julia
@model function fit_model(s₁, datasets)
    σ ~ InverseGamma(1, 2)
    A ~ Normal(sol.u.A, sol.u.A*0.1)
    a ~ Normal(sol.u.a, sol.u.a*0.1)
    b ~ Normal(sol.u.b, sol.u.b*0.1)
    ŝ₁ = map(data->ŝ(ψ, data, p_ExpLn(A, a, b)), datasets)
    ŝ₁ = vcat(ŝ₁...)
    for (index, ŝ) in enumerate(ŝ₁)
        s₁[index] ~ MvNormal([ŝ], σ^2*I)
    end
    return nothing
end
```

Next, we create a vector of the Uniaxial stresses from the datasets:

{cell}
```julia
s₁ = map(data->getindex.(data.data.s,1), kawabata_data)
s₁ = vcat(s₁...)
s₁ = map(s->[s], s₁)
```

Then we can create the model for fitting with `Turing.jl`:

{cell outputshow=false resultshow=false}
```julia
model = fit_model(s₁, kawabata_data);
```

Finally, sampling the model gives:

{cell output=false resultshow = false}
```julia
chns = sample(model, NUTS(0.65), MCMCSerial(), 3000, 2)
```

Finally, plotting the distributions

{cell outputshow = false}
```julia
params = names(chns, :parameters)
n_chains = length(chains(chns))
n_samples = length(chns)
f = Figure()
for (i, param) in enumerate(params)
    ax = Axis(f[i, 1], xlabel = string(param))
    colors = Makie.current_default_theme().attributes[:palette][][:color][]
    for chain in 1:n_chains
        values = chns[:, param, chain]
        density!(ax, values; color = (colors[chain], 0.3), strokecolor = :black, strokewidth = 1)
    end
end
f
```

Now, using the fitted distributions a data retrodiction plot can be generated to see how distributions impact the performance of the model:

{cell}
```julia
f = Figure()
ax = Axis(f[1,1], xlabel = "Stretch [-]", ylabel = "Stress [MPa]")
posterior_samples = sample(chns[[:A, :a, :b]], 2000, replace = false)
test_data = kawabata_data[2]
for posterior_p in eachrow(Array(posterior_samples))
    ŝ₂ = getindex.(Hyperelastics.predict(ψ, test_data, p_ExpLn(posterior_p[1], posterior_p[2],posterior_p[3])).data.s, 1)
    lines!(ax, getindex.(test_data.data.λ, 2), ŝ₂, alpha = 0.001, color = "#CCCCCC", strokewidth = 0.2)
end
scatter!(ax, getindex.(test_data.data.λ, 2), getindex.(test_data.data.s, 1), color = :black)
f
```