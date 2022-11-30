using Hyperelastics
using AbstractDifferentiation, ForwardDiff, FiniteDifferences
using Optimization, OptimizationOptimJL, ComponentArrays
using Plots
using QuadGK, Integrals, FastGaussQuadrature, HCubature
using Interpolations, DataInterpolations
using LinearAlgebra, Tullio

## Test Stretches
λ₁ = [1.04, 2.2, 3.1, 3.4]
data = Kawabata1981.(λ₁)
λ⃗ = vcat(map(d -> collect(d.λ⃗), data)...)
λ₂ = getindex.(collect.(λ⃗), 2)
s⃗ = vcat(map(d -> collect(d.s⃗), data)...)
data = biaxial_data(getindex.(s⃗, 1), getindex.(s⃗, 2), getindex.(λ⃗, 1), getindex.(λ⃗, 2))

n1 = 15
p₀ = rand(n1)
function Pchain(p)
    itp = interpolate(p, BSpline(Cubic(Periodic(OnGrid()))))
    sitp = scale(itp, range(1 / maximum(λ₂)^2, maximum(λ₂), length=n1))
    return sitp
end

# Set the model Parameters
ps = (
    Pchain=Pchain,
    p₀=p₀,
    weights=ones(length(data.s⃗)),
    coeff_loss=x -> 0.0,
    solver=GradientDescent(),
    data=data
)

# Create the Model
# pchain, λchains, s̃ᵢ, ∂Ψ∂λᵢ, ∂Ψ̂∂λᵢ = MacroMicroMacro(ps)
W = MacroMicroMacro(ps);

## Compare the Experimental and Predicted Stresses
s₂ = s̃ᵢ.(λ⃗, 2)
s₁ = s̃ᵢ.(λ⃗, 1)
σ₁ = getindex.(data.s⃗, 1)
σ₂ = getindex.(data.s⃗, 2)
scatter(getindex.(λ⃗, 2), σ₂)
scatter!(getindex.(λ⃗, 2), s₂)

scatter(getindex.(λ⃗, 2), σ₁)
scatter!(getindex.(λ⃗, 2), s₁)
## Test Set
λ₁ = 1.3
data = Kawabata1981(λ₁)
λ₂_test = getindex.(data.λ⃗, 2)
λ⃗ = map(x -> [x[1], x[2], 1 / x[1] / x[2]], zip(fill(λ₁, 30), range(extrema(λ₂_test)..., length=30)))
s₁_test = s̃ᵢ.(λ⃗, 1)
s₂_test = s̃ᵢ.(λ⃗, 2)
σ₁_test = getindex.(data.s⃗, 1)
σ₂_test = getindex.(data.s⃗, 2)
plot(getindex.(λ⃗, 2), s₂_test)
scatter!(λ₂_test, σ₂_test, legend=false)

plot!(getindex.(λ⃗, 2), s₁_test)
scatter!(λ₂_test, σ₁_test, legend=false)

## Plot the Stress-Stretch Relationship Determined for the Chain
plot(
    range(extrema(pchain.ranges[1])..., length=100),
    x -> pchain(x) / 1e6,
    xlabel="λ-chain",
    ylabel="P-chain [MPa]",
    legend=false
)
scatter!(
    pchain.ranges[1],
    pchain.itp.coefs .* 1e-6
)
