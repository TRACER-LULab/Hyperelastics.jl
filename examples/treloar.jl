# # Package Imports
using Hyperelastics
using IntervalArithmetic, IntervalOptimisation, LossFunctions
using Optimization, OptimizationOptimJL, ComponentArrays, DataInterpolations
using CairoMakie, ColorSchemes
using MKL
# # Treloar's UnPiaxial Data
s₁ = [0.0, 0.2856, 0.3833, 0.4658, 0.5935, 0.6609, 0.8409, 1.006, 1.2087, 1.5617, 1.915, 2.2985, 2.6519, 3.0205, 3.3816, 3.7351, 4.0812, 4.4501, 4.8414, 5.2026, 5.5639] * 1e6
λ₁ = [1.0, 1.4273, 1.6163, 1.882, 2.1596, 2.4383, 3.0585, 3.6153, 4.1206, 4.852, 5.4053, 5.7925, 6.1803, 6.4787, 6.6627, 6.936, 7.133, 7.1769, 7.2712, 7.4425, 7.512]
λ⃗_predict = collect(map(x -> [x, 1 / √x, 1 / √x], range(minimum(λ₁), maximum(λ₁), length=30)))
# # Create a Uniaxaial Test Results Object
data = UniaxialHyperelasticData(s₁, λ₁)

# ## Fit the Gent Model
# $W(\vec{\lambda}) = -\frac{\mu J_m}{2}\log{\bigg(1-\frac{I_1-3}{J_m}\bigg)}$
#
# Initial guess for the parameters
ψ = Gent()
p₀ = ComponentVector(μ=240e3, Jₘ=30.0)
# Create the optimization problem and solve
HEProblem = HyperelasticProblem(
    data,
    ψ,
    p₀,
    [],
)
sol = solve(HEProblem, LBFGS())
# $\mu$ = 240kPa, $J_m$ = 79.97
# Predict the new stresses
function ŝ(λ⃗, ψ, p)
    s⃗ = NominalStressFunction(ψ, λ⃗, p)
    return s⃗ .- s⃗[3] .* λ⃗[3] ./ λ⃗
end
ŝ₁ = map(λ -> ŝ(λ, ψ, sol.u)[1], λ⃗_predict)
# Plot the Results
fig = Figure(font = "CMU Serif")
ax = Makie.Axis(
    fig[1, 1],
    xlabel="Stretch",
    ylabel="Stress [MPa]",
    palette=(color=ColorSchemes.Egypt,),
    )
scatter!(
    ax,
    getindex.(data.λ⃗, 1),
    getindex.(data.s⃗, 1) ./ 1e6,
    marker='◆',
    markersize=25,
    color=:black,
    label="Experimental"
)

lines!(
    ax,
    getindex.(λ⃗_predict, 1),
    ŝ₁ ./ 1e6,
    label="Predicted $(string(ψ)[1:end-2])"
)
current_figure()
save("examples/" * string(ψ)[1:end-2] * ".png", current_figure()) #src
# ![Gent Plot](../examples/gent.png)

# ## Using the EdwardVilgis Model
ψ = EdwardVilgis()
p₀ = ComponentVector(Ns=20e3, Nc=20e3, α=0.05, η=0.05)
HEProblem = HyperelasticProblem(
    data,
    ψ,
    p₀,
    []
)
sol = solve(HEProblem, LBFGS())
#
# Plot and compare the stresses
ŝ₁ = map(λ -> ŝ(λ, ψ, sol.u)[1], λ⃗_predict)
lines!(
    ax,
    getindex.(λ⃗_predict, 1),
    ŝ₁ ./ 1e6,
    label="Predicted $(string(ψ))"
)
current_figure()
# ## Using the ModifiedFloryErman Model
ψ = ModifiedFloryErman()
p₀ = ComponentVector(μ=240e3, N=60.0, κ=100.0)
HEProblem = HyperelasticProblem(
    data,
    ψ,
    p₀,
    []
)
sol = solve(HEProblem, LBFGS())
#
# Plot and compare the stresses
ŝ₁ = map(λ -> ŝ(λ, ψ, sol.u)[1], λ⃗_predict)
lines!(
    getindex.(λ⃗_predict, 1),
    ŝ₁ ./ 1e6,
    label="Predicted ModifiedFloryErman"
)
current_figure()
# ## Using the MCC Model
ψ = MCC()
p₀ = ComponentVector(κ=100000000.0, μkT=10e3, ζkT=10e3)
HEProblem = HyperelasticProblem(
    data,
    ψ,
    p₀,
    []
)
sol = solve(HEProblem, LBFGS())
#
# Plot and compare the stresses
ŝ₁ = map(λ -> ŝ(λ, ψ, sol.u)[1], λ⃗_predict)
lines!(
    ax,
    getindex.(λ⃗_predict, 1),
    ŝ₁ ./ 1e6,
    label="Predicted $(string(ψ))"
)
current_figure()
# ## Using the NeoHookean Model
# $W(\vec{\lambda}) = \frac{\mu}{2}(I_1-3)$
ψ = NeoHookean()
p₀ = ComponentVector(μ=200e3)
HEProblem = HyperelasticProblem(
    data,
    NeoHookean(),
    p₀,
    []
)
sol = solve(HEProblem, LBFGS())
# $\mu$ = 534kPa
#
# Plot and compare the stresses
ŝ₁ = map(λ -> ŝ(λ, ψ, sol.u)[1], λ⃗_predict)
lines!(
    ax,
    getindex.(λ⃗_predict, 1),
    ŝ₁ ./ 1e6,
    label="Predicted NeoHookean"
)
save("examples/" * string(ψ)[1:end-2] * ".png", current_figure()) #src
# ![Neohookean Plot](../examples/NeoHookean.png)
# ## Sussman-Bathe Model
# $W(\vec{\lambda}) = \sum\limits_{i=1}^{3} w(\lambda_i)$
ψ = SussmanBathe(data, 5, DataInterpolations.QuadraticSpline)
ŝ₁ = map(λ -> ŝ(λ, ψ, sol.u)[1], λ⃗_predict)

lines!(
    ax,
    getindex.(λ⃗_predict, 1),
    ŝ₁ ./ 1e6,
    label="Predicted NeoHookean"
)

current_figure()

savefig("examples/sussmanbathe.png") #src
# ![Sussman Bathe Plot](../examples/sussmanbathe.png)
plot!() #src
# # Using Turing.jl for Parameter Estimation
using Turing, Distributions, LinearAlgebra
using MCMCChains, SciMLExpectations, KernelDensity
using Cuba
Turing.setadbackend(:forwarddiff)
function Makie.plot(ch::Chains)
    fig = Figure()
    for (ind, param) in enumerate(ch.name_map.parameters)
        ax = Makie.Axis(fig[ind, 1], xlabel=string(param))
        for (ind2, datavec) in enumerate(eachcol(getindex(ch, param).data))
            # Get current default colorpalette
            colors = Makie.current_default_theme().attributes[:palette][][:color][]
            Makie.density!(ax, datavec, color=(:black, 0.0),
                strokearound=true,
                strokewidth=2,
                strokecolor=colors[ind2%length(colors)]
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
Jₘ_min = maximum(I₁.(collect.(data.λ⃗))) - 3
function ŝ(λ⃗, ψ, p)
    s⃗ = NominalStressFunction(ψ, λ⃗, p)
    return s⃗ .- s⃗[3] .* λ⃗[3] ./ λ⃗
end
@model function fitHE(s₁, data)
    ## Prior Distributions
    σ ~ InverseGamma(1, 2) # noise in the measurement data
    μ ~ truncated(Normal(2.4e5, 0.5e5), 2.2e5, 2.6e5) # Normal for μ
    Jₘ ~ truncated(Normal(300.0, 30.0), lower=Jₘ_min) # Truncated Normal for Jₘ with lower bound

    ## Simulate the data
    ŝ₁ = map(λ -> ŝ(λ, Gent(), ps(μ, Jₘ))[1], data.λ⃗)

    ## Observations
    for (index, ŝ) in enumerate(ŝ₁)
        s₁[index] ~ MvNormal([ŝ], σ^2 * I)
    end

    return nothing
end
test_s = map(s -> [s], s₁)
model = fitHE(test_s, data)
# # Samble the distributions to fit the data and print the results
chain = sample(model, NUTS(0.65), MCMCThreads(), 6000, 4, progress=true)
# $\mu$ = 245kPa ± 5.238kPa, $J_m$ = 80.9±1.1583
Makie.plot(chain)
save("examples/chain.png", current_figure()) #src
# ![chain](../examples/chain.png)
# Data Retrodiction to observe the results with 300 samples from the chain
fig = Figure()
ax = Makie.Axis(fig[1, 1], xlabel="Stretch", ylabel="Stress [MPa]")
posterior_samples = sample(chain[[:μ, :Jₘ]], 1000; replace=false)
for p in eachrow(Array(posterior_samples))
    ŝ₁ = map(λ -> ŝ(λ, Gent(), (μ=p[1], Jₘ=p[2]))[1], λ⃗_predict)
    lines!(ax, getindex.(λ⃗_predict, 1), ŝ₁ ./ 1e6, alpha=0.1, color="#BBBBBB")
end
scatter!(ax, λ₁, s₁ ./ 1e6, color=:black)
current_figure()
save("examples/dataretrodiction.png", current_figure()) #src
# ![retrodiction](../examples/dataretrodiction.png)
# # Analysing the final distribution at the last data point
fig = Figure()
ax = Makie.Axis(fig[1, 1], xlabel="Stress [MPa]", ylabel="Probability")
posterior_samples = sample(chain[[:μ, :Jₘ]], 2000; replace=false)
μ_dist = kde(vec(Array(chain[:μ])))
Jₘ_dist = kde(vec(Array(chain[:Jₘ])))
s(x) = ŝ(λ⃗_predict[end], Gent(), ps(μ_dist(1), x[2]))[1]
s₁_values = map(x-> s(x), eachrow(Array(posterior_samples)))
density!(ax, s₁_values, color=(:slategray, 0.4))
current_figure()
# # Generating partial derivatives of the SEF with Symbolics.jl
using Symbolics
# ## Create the symbolic variables required
@syms μ Jₘ
@syms λ₁ λ₂ λ₃
# ## Make the symbolic model
W = StrainEnergyDensityFunction(Gent(), [λ₁, λ₂, λ₃], (μ=μ, Jₘ=Jₘ))
# ## Create the partial derivative operators for the principal stretches
∂λ = Differential.([λ₁, λ₂, λ₃])
# ## Differentiate the SEF with respect to the principal stretches
∂W∂λ = map(∂ -> ∂(W), ∂λ)
∂W∂λ = expand_derivatives.(∂W∂λ)
# ## Create a function from the symbolic expression for each partial derivative
# Using the mean values from the Bayesian Parameter Estimation as the material properties
∂W∂λ = substitute.(∂W∂λ, (Dict(μ => mean(chain[:μ]), Jₘ => mean(chain[:Jₘ])),))
∂W∂λ = simplify.(∂W∂λ)
∂W∂λ = map(∂ -> build_function(∂, [λ₁, λ₂, λ₃], expression=Val{false}), ∂W∂λ)
# Predict the results
s⃗_predict = map(λ -> map(∂ -> ∂(λ), ∂W∂λ), λ⃗_predict)
σ⃗_predict = map(x -> x[1] .* x[2], zip(s⃗_predict, λ⃗_predict))
σ₁ = getindex.(σ⃗_predict, 1) - getindex.(σ⃗_predict, 3)
s1 = σ₁ ./ getindex.(λ⃗_predict, 1)
# Plot the results
lines!(ax, getindex.(λ⃗_predict, 1), s1 ./ 1e6, color=:red, label="Symbolic mean")
current_figure()
save("examples/symbolic_plot.png", current_figure()) #src
# ![symbolic](../examples/symbolic_plot.png)
