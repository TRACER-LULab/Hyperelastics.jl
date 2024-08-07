# # Package Imports
using Hyperelastics
# using Unitful
using DifferentiationInterface
import ForwardDiff, FiniteDiff
using Optimization
using OptimizationOptimJL
using ComponentArrays
using CairoMakie
# # Load the Treloar 1994 Uniaxial Dataset
# treloar_data = [Treloar1944Uniaxial(),Treloar1944Uniaxial()]
treloar_data = Treloar1944Uniaxial()
λ₁ = getindex.(treloar_data.data.λ, 1)
s₁ = getindex.(treloar_data.data.s, 1)
treloar_data = HyperelasticUniaxialTest(λ₁, s₁, name = "test")

# ## Fit the Gent Model
# $W(\vec{\lambda}) = -\frac{\mu J_m}{2}\log{\bigg(1-\frac{I_1-3}{J_m}\bigg)}$
#
## Initial guess for the parameters
models = Dict(
    Gent => ComponentVector((μ = 240e-3, J_m = 79.0)),
    # EdwardVilgis => ComponentVector(Ns=0.10, Nc=0.20, α=0.001, η=0.001),
    # ModifiedFloryErman => ComponentVector(μ=0.24, N=50.0, κ=10.0),
    # NeoHookean => ComponentVector(μ=200e-3),
    # NonaffineMicroSphere => HyperelasticsZygoteExComponentVector(μ=0.292, N=22.5, p=1.471, U=0.744, q=0.1086),
    # Beda => ComponentVector(C1=0.1237, C2=0.0424, C3=7.84e-5, K1=0.0168, α=0.9, β=0.68, ζ=3.015)
)

#
# f = Figure()
# ax = Makie.Axis(f[1, 1], xlabel="Stretch", ylabel="Stress [MPa]")

# scatter!(ax, getindex.(pred.data.λ, 1), getindex.(treloar_data.data.s, 1), label="Treloar Data")
for (ψ, p₀) in models
    prob = HyperelasticProblem(ψ(), treloar_data, p₀, ad_type = AutoForwardDiff())
    display(prob.u0)
    sol = solve(prob, NelderMead())
    @show sol
    pred = predict(ψ(), treloar_data, p₀, ad_type = AutoForwardDiff())
    # @show pred
    # lines!(ax, getindex.(pred.data.λ, 1), getindex.(pred.data.s, 1), label=string(ψ))
end
# axislegend(position=:lt)
# f
## # For multiple tests
f = Figure()
ax = Makie.Axis(f[1, 1], xlabel = "Stretch", ylabel = "Stress [MPa]")
kawabata_data = map(
    λ₁ -> Kawabata1981(λ₁),
    [
        1.040,
        1.060,
        1.080,
        1.100,
        1.120,
        1.14,
        1.16,
        1.2,
        1.24,
        1.3,
        1.6,
        1.9,
        2.2,
        2.5,
        2.8,
        3.1,
        3.4,
        3.7,
    ],
)
scatter!(
    ax,
    getindex.(kawabata_data[1].data.λ, 2),
    getindex.(kawabata_data[1].data.s, 1),
    label = "Treloar Data",
)
for (ψ, p₀) in models
    HEProblem = HyperelasticProblem(ψ(), kawabata_data, p₀)
    sol = solve(HEProblem, LBFGS())
    pred = predict(ψ(), kawabata_data, sol.u)
    lines!(
        ax,
        getindex.(pred[1].data.λ, 2),
        getindex.(pred[1].data.s, 1),
        label = string(ψ),
    )
end
f
save("model_examples.png", f)
# ## Sussman-Bathe Model
using DataInterpolations
# $W(\vec{\lambda}) = \sum\limits_{i=1}^{3} w(\lambda_i)$
ψ = SussmanBathe(treloar_data, 5, DataInterpolations.QuadraticSpline)
pred = predict(ψ, treloar_data, [])
ŝ₁ = getindex.(pred.data.s, 1)
lines!(ax, λ₁, ŝ₁, label = "Sussmanbathe")
current_figure()
axislegend(position = :lt)
f
savefig("examples/sussmanbathe.png") #src
# ![Sussman Bathe Plot](../examples/sussmanbathe.png)
plot!() #src
# # Using Turing.jl for Parameter Estimation
using Turing, Distributions, LinearAlgebra
using MCMCChains, SciMLExpectations, KernelDensity
Turing.setadbackend(:forwarddiff)
function Makie.plot(ch::Chains)
    fig = Figure()
    for (ind, param) in enumerate(ch.name_map.parameters)
        ax = Makie.Axis(fig[ind, 1], xlabel = string(param))
        for (ind2, datavec) in enumerate(eachcol(getindex(ch, param).data))
            # Get current default colorpalette
            colors = Makie.current_default_theme().attributes[:palette][][:color][]
            Makie.density!(
                ax,
                datavec,
                color = (:black, 0.0),
                strokearound = true,
                strokecolor = colors[ind2%length(colors)],
            )
        end
    end
    display(fig)
    return fig
end
# Create a type to store the parameters used to simulate the data.
struct ps{T}
    μ::T
    Jₘ::T
end
Jₘ_min = maximum(I₁.(data.data.λ)) - 3
function ŝ(ψ, test, p)
    ŷ = predict(ψ, test, p)
    s = getindex.(ŷ.data.s, 1)
    return s
end
@model function fitHE(s₁, data)
    ## Prior Distributions
    σ ~ InverseGamma(1, 2) # noise in the measurement data
    μ ~ truncated(Normal(0.24, 0.05), 0.2, 0.3) # Normal for μ
    Jₘ ~ truncated(Normal(80.0, 10.0), lower = Jₘ_min) # Truncated Normal for Jₘ with lower bound

    ## Simulate the data
    ŝ₁ = ŝ(Gent(), data, ps(μ, Jₘ))
    # ŝ₁ = map(λ -> ŝ(λ, Gent(), ps(μ, Jₘ))[1], data.λ⃗)

    ## Observations
    for (index, ŝ) in enumerate(ŝ₁)
        s₁[index] ~ MvNormal([ŝ], σ^2 * I)
    end

    return nothing
end
test_s = map(s -> [s], s₁)
model = fitHE(test_s, data)
# # Samble the distributions to fit the data and print the results
chain = sample(model, NUTS(0.65), MCMCThreads(), 6000, 4, progress = true)
# $\mu$ = 245kPa ± 5.238kPa, $J_m$ = 80.9±1.1583
Makie.plot(chain)
save("examples/chain.png", current_figure()) #src
# ![chain](../examples/chain.png)
# Data Retrodiction to observe the results with 300 samples from the chain
fig = Figure()
ax = Makie.Axis(fig[1, 1], xlabel = "Stretch", ylabel = "Stress [MPa]")
posterior_samples = sample(chain[[:μ, :Jₘ]], 1000; replace = false)
for p in eachrow(Array(posterior_samples))
    ŝ₁ = map(λ -> ŝ(λ, Gent(), (μ = p[1], Jₘ = p[2]))[1], λ⃗_predict)
    lines!(ax, getindex.(λ⃗_predict, 1), ŝ₁ ./ 1e6, alpha = 0.1, color = "#BBBBBB")
end
scatter!(ax, λ₁, s₁ ./ 1e6, color = :black)
current_figure()
save("examples/dataretrodiction.png", current_figure()) #src
# ![retrodiction](../examples/dataretrodiction.png)
# # Analysing the final distribution at the last data point
fig = Figure()
ax = Makie.Axis(fig[1, 1], xlabel = "Stress [MPa]", ylabel = "Probability")
posterior_samples = sample(chain[[:μ, :Jₘ]], 2000; replace = false)
μ_dist = kde(vec(Array(chain[:μ])))
Jₘ_dist = kde(vec(Array(chain[:Jₘ])))
s(x) = ŝ(λ⃗_predict[end], Gent(), ps(μ_dist(1), x[2]))[1]
s₁_values = map(x -> s(x), eachrow(Array(posterior_samples)))
density!(ax, s₁_values, color = (:slategray, 0.4))
current_figure()
