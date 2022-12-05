using Hyperelastics
using Optimization, OptimizationOptimJL, ComponentArrays
using CairoMakie
using Interpolations, DataInterpolations
using LinearAlgebra, Tullio

## Test Stretches
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
