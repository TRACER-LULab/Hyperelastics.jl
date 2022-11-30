# # Package Imports
using Hyperelastics
using Optimization, OptimizationOptimJL, ComponentArrays
using Plots
using AbstractDifferentiation
using QuadGK, DataInterpolations
using LinearAlgebra
using NonlinearSolve
pgfplotsx() # src
## # Wavy Data from Sussman and Bathe equation (11)
e = range(-2, 2, length=100)
λ = exp.(e)
τ = @. (exp(3e) - exp(-3e)) * (1 + 0.2sin(10e)) - (exp(-1.5e) - exp(1.5e)) * (1 + 0.2sin(-5e))
τ_noise = τ .+ randn(100) * 30
s = τ ./ λ
data = uniaxial_data(s, λ)

## Curvature Contrained Susman Bathe
σ = Hyperelastics.StabilitySmoothedSussmanBathe((s⃗=s, λ⃗=λ, k=3, n=50, q=0.0))
e⃗ = range(-2, 2, length=100)
σ⃗ = σ.(e⃗)
plot(e⃗, σ⃗)
scatter!(e, τ_noise, markersize=2)
##
# # B-Spline Interpolation
function N(t; k, i, t⃗)
    if k ≥ 2
        return (t - t⃗[i]) / (t⃗[i+k-1] - t⃗[i]) * N(t, k=k - 1, i=i, t⃗=t⃗) + (t⃗[i+k] - t) / (t⃗[i+k] - t⃗[i+1]) * N(t, k=k - 1, i=i + 1, t⃗=t⃗)
    elseif t⃗[i] ≤ t < t⃗[i+1]
        return 1.0
    else
        return 0.0
    end
end
# # Regression Fitting
f̂(B, p) = [1 / 2 / p.N * (p.N̂ * B - τ)' * I * (p.N̂ * B - τ)]
# # Curvature Based Smoothing
function f̃_curve(B, p)

end
#
k = 3
n = 101
t⃗ = range(extrema(e)..., length=n + k + 1)
N⃗(E) = [N(E, k=k, i=i, t⃗=t⃗) for i in 1:n]'
N̂ = vcat(N⃗.(e)...)
func = OptimizationFunction(f̂, Optimization.AutoForwardDiff())
prob = OptimizationProblem(func, randn(n), (N=n, N̂=N̂))
sol = solve(prob, LBFGS())
#######
B = sol.u
Ê = e
𝐍(E) = [N(E, k=k, i=i, t⃗=t⃗) for i in 1:n]'
# E = [zeros(k);Ê;zeros(k)]
e⃗ = range(-2, 2, length=100)
sigma = vcat(𝐍.(e⃗)...) * B
plot(e⃗, sigma, markersize=2)
scatter!(e, τ, markersize=2)
#######

# # Sussman-Bathe Model
# ```math
# W(\vec{\lambda}) = \sum\limits_{i = 1}^3w(\lambda_i)
# ```
# where $w^\prime(\lambda) = \sum\limits_{j = 0}^{k}s_1\Big(\lambda^{(4^{-i})}\Big)+s_1\Big(\lambda^{(-\frac{1}{2}4^{-i})}\Big)$ such that $s_1(\lambda)$ is an interpolation of the stress-stretch data provided and $k$ is the limit of the summation. In the original paper, $k=4$ or $5$ was sufficient for accurate results.
+W = Hyperelastics.SussmanBathe((s⃗=data.s⃗, λ⃗=data.λ⃗, k=3))
adb = AD.FiniteDifferencesBackend()
ŝ = s⃗̂(W, collect.(data.λ⃗), adb=adb)
ŝ₁ = getindex.(ŝ, 1)
plot(λ, (ŝ₁ .* λ - s .* λ) ./ (s .* λ) .* 100, xaxis="λ", yaxis="Relative Error")
## Smoothed Sussman Bathe
Ws = StabilitySmoothedSussmanBathe((s⃗=data.s⃗, λ⃗=data.λ⃗, k=3, n=50))

##
We = Hyperelastics.LatoreeMontans(
    (
        type=:T1,
        σ̃₁=τ,
        Ẽ₁=e,
        Ẽ₂=-e / 2,
        k=3
    ),
    (
        type=:none,
    )
)

adb = AD.FiniteDifferencesBackend();
E = zip(e, -e / 2, -e / 2)
E = collect.(E)
E = map(x -> [x; [0.0, 0.0, 0.0, 0.0]], E)
∂We(x) = AD.gradient(adb, We, x)[1]
s⃗ = ∂We.(E)
s₁ = getindex.(s⃗, 1) - getindex.(s⃗, 3)
# ŝ = s⃗̂(We, L, adb=adb)
# ŝ₁ = getindex.(ŝ, 1)
plot(λ, s .* λ, xscale=:ln)
plot!(λ, ŝ₁ .* λ)
# # Test with data from Diani et al.
λ₁ = λ₃ = [1.00000, 1.00197, 1.09500, 1.09895, 1.19585, 1.19779, 1.29667, 1.29858, 1.39747, 1.39934, 1.49826, 1.50209, 1.59708, 1.59891, 1.69786, 1.69968, 1.79845, 1.79866, 1.89748, 1.89922, 1.99829, 1.99999, 2.09912, 2.10472, 2.19797, 2.20155, 2.29881, 2.30036, 2.39969, 2.40115]
e₁ = e₃ = log.(λ₁)

s₁_cal = [0.0000, 0.0129, 0.5653, 0.5840, 0.9316, 0.9372, 1.1986, 1.2034, 1.4456, 1.4500, 1.6728, 1.6815, 1.8901, 1.8937, 2.0776, 2.0813, 2.3142, 2.3147, 2.5519, 2.5561, 2.8188, 2.8238, 3.1255, 3.1424, 3.4224, 3.4337, 3.7788, 3.7850, 4.2247, 4.2314]
σ₁_cal = s₁_cal .* λ₁

s₃_trans = [0.0000, 0.0129, 0.5382, 0.5553, 0.8280, 0.8322, 1.0361, 1.0395, 1.1941, 1.1971, 1.3679, 1.3746, 1.5294, 1.5323, 1.6776, 1.6799, 1.7978, 1.7981, 1.9625, 1.9654, 2.1202, 2.1230, 2.2999, 2.3104, 2.4905, 2.4978, 2.7020, 2.7052, 2.8998, 2.9026]
σ₃_trans = s₃_trans .* λ₁

We = Hyperelastics.LatoreeMontans(
    (
        type=:T2B,
        σ̃₁=σ₁_cal,
        σ̃₃=σ₃_trans,
        Ẽ₁=e₁,
        Ẽ₃=e₃,
        Ê₂=(x, p) -> p[1] * x^3 + p[2] * x^2 + p[3] * x,
        p₀=[1.0, 1.0, 1.0],
        solver=LBFGS(),
        k=3
    ),
    (
        type=:none,
    )
)
#######
diff = ∂Ψ∂λᵢ.(λ⃗) .- ∂Ψ̂∂λᵢ.(λ⃗)
plot(getindex.(λ⃗, 2), getindex.(diff, 1))
plot!(getindex.(λ⃗, 2), getindex.(diff, 2))
plot!(getindex.(λ⃗, 2), getindex.(diff, 3))

##
x, w = gausslobatto(30)
W1(λ⃗) = sum(i->w[i] .* ∂Ψ∂λᵢ([λ⃗[1] * (x[i] + 1) / 2 + (1 - x[i]) / 2, λ⃗[2], λ⃗[3]])[1], eachindex(w))
W2(λ⃗) = sum(i->w[i] .* ∂Ψ∂λᵢ([λ⃗[1], λ⃗[2] * (x[i] + 1) / 2 + (1 - x[i]) / 2, λ⃗[3]])[2], eachindex(w))
W3(λ⃗) = sum(i->w[i] .* ∂Ψ∂λᵢ([λ⃗[1], λ⃗[2], λ⃗[3] * (x[i] + 1) / 2 + (1 - x[i]) / 2])[3], eachindex(w))

Δ(λ⃗) = [W1(λ⃗), W2(λ⃗), W3(λ⃗)]
as = Δ.(λ⃗)
plot(getindex.(as, 1), getindex.(as, 2), getindex.(as, 3), size = (500, 500))

∂W̃∂λ(λ⃗) = AD.gradient(AD.ForwardDiffBackend(), W, λ⃗)[1]
∂W̃∂λ([2.0, 1 / sqrt(2), 1 / sqrt(2)])
s̃ᵢ_W(λ⃗, i) = ∂W̃∂λ(λ⃗)[i] - ∂W̃∂λ(λ⃗)[3] * λ⃗[3] / λ⃗[i]
##
chain_plot = plot(range(extrema(λ₂)..., length=1000), pchain, xlims=(0, 3.2), xticks=0:0.8:3.2)

s₂ = Base.Fix2(s̃ᵢ_W, 2).(λ⃗)
s̃₂ = Base.Fix2(s̃ᵢ, 2).(λ⃗)
scatter(λ₂, σ₂)
plot(getindex.(λ⃗, 2), s₂)
plot!(getindex.(λ⃗, 2), s̃₂)

λ₂_test = range(extrema(λ₂)..., length=20)
plot(λ₂_test, s̃.(3.1, λ₂_test, 2))
plot(λ₂_test, s̃.(3.1, λ₂_test, 1))
scatter!(λ₂, σ₁)

n1 = 20
xs = collect(range(0.5, 3, length=n1))
ys = collect(range(0.5, 3, length=n1))
x_grid = [x for x = xs for y = ys]
y_grid = [y for x = xs for y = ys]
zs = s̃.(xs .* ones(n1)', ys' .* ones(n1), 2) ./ 1e6
surface(xs, ys, zs, xlabel="l1", ylabel="l2", zlims=(0.6, 1.3), size=(900, 900), camera=(45, 15))
plot!(λ₂, σ₂)

s₁ = p1(λ⃗)
s₂ = p2(λ⃗)
stress_plot = scatter(λ₂, σ₂)
plot!(getindex.(λ⃗, 2), s₂)
scatter(λ₂, σ₁)
plot!(getindex.(λ⃗, 2), s₁)

using AbstractDifferentiation
using ForwardDiff
using FiniteDifferences
s⃗(λ⃗) = AD.gradient(AD.FiniteDifferencesBackend(), Ψ, λ⃗)
s̃ = s⃗.(λ⃗)
s̃ = getindex.(s̃, 1)
s̃₃ = getindex.(s̃, 3)
s̃₂ = getindex.(s̃, 2)
λ̃₂ = getindex.(λ⃗, 2)
λ̃₃ = getindex.(λ⃗, 3)
s2 = @. s̃₂ - λ̃₃ / λ̃₂ * s̃₃

plot(λ̃₂, s2)
plot!(λ̃₂, s̃₂)
plot!(λ₂, σ₂)
