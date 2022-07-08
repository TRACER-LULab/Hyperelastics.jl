using Hyperelastics
using FastGaussQuadrature
using AbstractDifferentiation, ForwardDiff, FiniteDifferences
using DataInterpolations, Interpolations
using QuadGK
using Optimization, OptimizationOptimJL
using Plots
## --------------------
λ₁max = 6.0
λ₁ = range(1.0, λ₁max, length=101)
λ₂ = range(1.0, λ₁max, length=101)
λ⃗ = map(λ -> [λ[1], λ[2], 1 / λ[1] / λ[2]], zip(λ₁, λ₂))
WAB = ArrudaBoyce((μ=0.27e6, N=26.5))
s₁ = getindex.(s⃗̂(WAB, λ⃗), 1)
s₂ = getindex.(s⃗̂(WAB, λ⃗), 2)
plot(λ₁, s₁, markersize=2)
## -----
WAvgChain = DataDrivenAverageChainBehavior(
    (
    data=BiaxialHyperelasticData(s₁, s₂, λ₁, λ₁),
    pchain=(x, y) -> BSplineInterpolation(y, x, 3, :Uniform, :Uniform)
)
)

s1 = getindex.(s⃗̂(WAvgChain, λ⃗[5:end-5], adb=AD.FiniteDifferencesBackend()), 1)
plot!(λ₁[5:end-5], s1)
