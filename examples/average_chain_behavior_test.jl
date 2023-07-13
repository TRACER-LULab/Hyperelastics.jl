using Hyperelastics
using FastGaussQuadrature
using AbstractDifferentiation, ForwardDiff, FiniteDifferences
using DataInterpolations, Interpolations
using QuadGK
using Optimization, OptimizationOptimJL
using Plots
## --------------------
λ₁max = 6.0
λ₁ = range(1.0, λ₁max, length = 101)
λ₂ = range(1.0, λ₁max, length = 101)
λ⃗ = map(λ -> [λ[1], λ[2], 1 / λ[1] / λ[2]], zip(λ₁, λ₂))
WAB = ArrudaBoyce((μ = 0.27e6, N = 100.5))
s₁ = getindex.(s⃗̂(WAB, λ⃗), 1)
s₂ = getindex.(s⃗̂(WAB, λ⃗), 2)
plot(λ₁, s₁, markersize = 2)
## -----
struct parameters
    data::AbstractHyperelasticData
    pchain::Function
end
data = BiaxialHyperelasticData(s₁, s₂, λ₁, λ₂)
ps = parameters(data, (x, y) -> BSplineInterpolation(y, x, 4, :Uniform, :Uniform))
WAvgChain = DataDrivenAverageChainBehavior(ps)
## -----
λ₁_test = range(1.0, λ₁max, length = 301)
λ₂_test = range(1.0, λ₁max, length = 301)
λ⃗_test = map(λ -> [λ[1], λ[2], 1 / λ[1] / λ[2]], zip(λ₁_test, λ₂_test))
s1 = getindex.(s⃗̂(WAvgChain, λ⃗_test[10:end-10], adb = AD.FiniteDifferencesBackend()), 1)
plot!(λ₁_test[10:end-10], s1, legend = false)
## -----
using Surrogates
using SurrogatesRandomForest
bounds = [1.0, 1.0], [λ₁max, λ₁max]
f = x -> WAvgChain([x[1], x[2], 1 / x[1] / x[2]])
x_train = sample(6000, bounds..., SobolSample())
y_train = f.(x_train)
wend = RadialBasis(x_train, y_train, bounds[1], bounds[2])
y_test = Float64.(wend.(x_train))
scatter3d(getindex.(x_train, 1), getindex.(x_train, 2), y_train)
scatter3d!(getindex.(x_train, 1), getindex.(x_train, 2), y_test)
s⃗(x) = AD.gradient(AD.ForwardDiffBackend(), x -> wend([x[1], x[2]]), x)[1]
s1 = getindex.(s⃗.(λ⃗_test), 1)
plot!(λ₁_test, s1)
