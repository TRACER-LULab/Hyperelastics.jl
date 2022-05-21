module Hyperelastics
using Reexport
using LossFunctions
using GalacticOptim
using AbstractDifferentiation
using Tullio
using Reexport
using ComponentArrays
using SpecialFunctions

# Write your package code here.
export HyperelasticData, uniaxial_data, biaxial_data, HyperelasticProblem
export I₁, I₂, I₃, J, s⃗̂

include("BasicDefinitions.jl")
include("HyperelasticModels.jl")

export s⃗̂, I₁, I₂, I₃, J
struct HyperelasticData
    s⃗
    λ⃗
end

"""
biaxial_data(s₁, s₂, λ₁, λ₂)

Create a biaxial hyperelastic data object from arrays of test data. The function returns a HyperelasticData object with the stresses and principal stretches. Currently, this assumes the material is incompressible.
"""
function biaxial_data(s₁, s₂, λ₁, λ₂)
    s⃗ = zip(s₁, s₂)
    λ⃗ = zip(λ₁, λ₂, (λ₁ .* λ₂) .^ (-1))
    return HyperelasticData(s⃗, λ⃗)
end

"""
uniaxial_data(s₁, λ₁)

Create a uniaxial hyperelastic data object from arrays of test data. The function returns a HyperelasticData object with the stresses and principal stretches. Currently, this assumes the material is incompressible.
"""
function uniaxial_data(s₁, λ₁)
    s⃗ = zip(s₁)
    λ⃗ = zip(λ₁, (λ₁) .^ (-0.5), (λ₁) .^ (-0.5))
    return HyperelasticData(s⃗, λ⃗)
end

include("wypiwyg.jl")

"""
---
HyperelasticProblem(data::HyperelasticData, model, u₀, ps; loss=L2DistLoss(), agg=AggMode.Mean(), cons=(x, p) -> [true], kwargs...)
---
Returns an `OptimizationProblem` for solving with GalacticOptim.jl. `data` is the hyperelastic experimental data, `model` is the strain energy density as a function of the parameters (i.e. `f(p) = W(p)(λ⃗)`). `ps` is any hyperparameters for the model (currently not supported). `loss` defines the loss function to be used in the optimization. Currently defaults to the ``L^2``-norm between the predicted and experimental data. `agg` defines the aggregration mode of the errors, defaults to the mean of the errors. `cons` define any constrain equations involving the parameters of the model. `kwargs` are passed to `OptimizationProblem`. To set parameter bounds, use the keywords `lb` and `ub` respectively. 

"""
function HyperelasticProblem(data::HyperelasticData, model, u₀, ps; loss=L2DistLoss(), agg=AggMode.Mean(), cons=(x, p) -> [true], kwargs...)
    s = hcat(collect.(data.s⃗)...)

    stresses_provided = size(s, 1)

    s⃗(p) = s⃗̂(model(p), collect.(data.λ⃗))

    function ŝ(p)
        s₁₂₃ = s⃗(p)
        return hcat(s₁₂₃...)[1:stresses_provided, :]
    end

    f(p, _) = [value(loss, s, ŝ(p), agg)]
    func = OptimizationFunction(f, GalacticOptim.AutoForwardDiff(), cons=cons)
    return OptimizationProblem(func, u₀, ps; kwargs...)
end

"""
---
HyperelasticProblem(data::Vector{HyperelasticData}, model, u₀, ps; loss=L2DistLoss(), agg=AggMode.Mean(), cons=(x, p) -> [true], kwargs...)
---
Returns an `OptimizationProblem` for solving with GalacticOptim.jl. `data` is a vector of hyperelastic data and fits to all sets equally, `model` is the strain energy density as a function of the parameters (i.e. `f(p) = W(p)(λ⃗)`). `ps` is any hyperparameters for the model (currently not supported). `loss` defines the loss function to be used in the optimization. Currently defaults to the ``L^2``-norm between the predicted and experimental data. `agg` defines the aggregration mode of the errors, defaults to the mean of the errors. `cons` define any constrain equations involving the parameters of the model. `kwargs` are passed to `OptimizationProblem`. To set parameter bounds, use the keywords `lb` and `ub` respectively. 

"""
function HyperelasticProblem(data::Vector{HyperelasticData}, model, u₀, ps; loss=L2DistLoss(), agg=AggMode.Mean(), cons=(x, p) -> [true], kwargs...)

    s = getindex.(vcat(collect.(zip(getfield.(data, :s⃗)...))...), 1) |> transpose
    λ = collect.(vcat(collect.(zip(getfield.(data, :λ⃗)...))...))

    stresses_provided = size(s, 1)

    s⃗(p) = s⃗̂(model, p, λ)

    function ŝ(p)
        s₁₂₃ = s⃗(p)
        return hcat(s₁₂₃...)[1:stresses_provided, :]
    end

    f(p, _) = [value(loss, s, ŝ(p), agg)]
    func = OptimizationFunction(f, GalacticOptim.AutoForwardDiff(), cons=cons)
    return OptimizationProblem(func, u₀, ps; kwargs...)
end

end
