# # Package Imports 
using Hyperelastics
using GalacticOptim
using GalacticOptimJL
using ComponentArrays
using Plots
pgfplotsx()
# # Treloar's Uniaxial Data
s₁ = [0.0, 0.2856, 0.3833, 0.4658, 0.5935, 0.6609, 0.8409, 1.006, 1.2087, 1.5617, 1.915, 2.2985, 2.6519, 3.0205, 3.3816, 3.7351, 4.0812, 4.4501, 4.8414, 5.2026, 5.5639] * 1e6
λ₁ = [1.0, 1.4273, 1.6163, 1.882, 2.1596, 2.4383, 3.0585, 3.6153, 4.1206, 4.852, 5.4053, 5.7925, 6.1803, 6.4787, 6.6627, 6.936, 7.133, 7.1769, 7.2712, 7.4425, 7.512]
λ⃗_predict = collect(map(x -> [x, 1 / √x, 1 / √x], range(minimum(λ₁), maximum(λ₁), length=30)))
# # Create a Uniaxaial Test Results Object
data = uniaxial_data(s₁, λ₁)

# # Fit the Gent Model
# $W(\vec{\lambda}) = -\frac{\mu J_m}{2}\log{\bigg(1-\frac{I_1-3}{J_m}\bigg)}$
#
# Initial guess for the parameters
p₀ = ComponentVector(μ=1e5, Jₘ=55.0)
# Set the bounds for the value $J_m$
Jₘ_min = I₁(collect.(data.λ⃗)[end]) - 3
lb = ComponentVector(μ=1.0, Jₘ=Jₘ_min)
ub = ComponentVector(μ=Inf, Jₘ=Inf)
# Create the optimization problem and solve
HEProblem = HyperelasticProblem(
    data,
    Gent,
    p₀,
    [],
    lb=lb,
    ub=ub
)
sol = solve(HEProblem, LBFGS())
# $\mu$ = 240kPa, $J_m$ = 79.97
# 
# Predict the new stresses
W = Gent(sol.u)
ŝ = s⃗̂(W, λ⃗_predict)
ŝ₁ = getindex.(ŝ, 1)
# Plot the Results
scatter(
    getindex.(data.λ⃗, 1), 
    getindex.(data.s⃗, 1) ./ 1e6, 
    label="Experimental"
)
plot!(
    getindex.(λ⃗_predict,1), 
    ŝ₁ ./ 1e6, 
    label="Predicted Gent"
)
plot!(xlabel="Stretch", ylabel="Stress [MPa]", legend=:topleft) #src
savefig("examples/gent.png") #src
# ![Gent Plot](examples/gent.png)
# # Using the NeoHookean Model
# $W(\vec{\lambda}) = \frac{\mu}{2}(I_1-3)$
p₀ = ComponentVector(μ=100e3)
HEProblem = HyperelasticProblem(
    data,
    NeoHookean,
    p₀,
    []
)
sol = solve(HEProblem, LBFGS())
# $\mu$ = 534kPa
# 
# Plot and compare the stresses
W = NeoHookean(sol.u)
ŝ = s⃗̂(W, λ⃗_predict)
ŝ₁ = getindex.(ŝ, 1)
plot!(
    getindex.(λ⃗_predict, 1),
    ŝ₁ ./ 1e6,
    label="Predicted NeoHookean"
)
savefig("examples/neohookean.png") #src 
# ![Neohookean Plot](examples/neohookean.png)
# # Sussman-Bathe Model
# $W(\vec{\lambda}) = \sum\limits_{i=1}^{3} w(\lambda_i)$
# 
# Note: the Sussman-Bathe model currently only supports differentiation via FiniteDifferences.jl as the AbstractDifferentiation.jl backend
using FiniteDifferences
using AbstractDifferentiation

W = Hyperelastics.SussmanBathe((s⃗=data.s⃗, λ⃗=data.λ⃗, k=3))
ŝ = s⃗̂(
    W,
    λ⃗_predict,
    adb=AD.FiniteDifferencesBackend()
)

ŝ₁ = getindex.(ŝ, 1)

plot!(
    getindex.(λ⃗_predict, 1),
    ŝ₁ ./ 1e6,
    label="Predicted Sussman-Bathe k = 4"
)

savefig("examples/sussmanbathe.png") #src
# ![Sussman Bathe Plot](examples/sussmanbathe.png)
plot!() #src