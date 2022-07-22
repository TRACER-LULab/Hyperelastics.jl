module Hyperelastics
using Reexport
using LossFunctions
using Optimization
using AbstractDifferentiation, ForwardDiff
using Tullio
using Reexport
using SpecialFunctions

# Write your package code here.
export AbstractHyperelasticData, UniaxialHyperelasticData, BiaxialHyperelasticData, HyperelasticProblem
export I₁, I₂, I₃, J
export StrainEnergyDensityFunction, NominalStressFunction, TrueStressFunction

include("basic_definitions.jl")
include("datasets.jl")
include("data_driven.jl")
include("isotropic_incompressible_models.jl")
include("macro_micro_macro_model.jl")
include("average_chain_behavior.jl")

abstract type AbstractHyperelasticData end
abstract type AbstractHyperelasticModel end
"""
uniaxial_data(s₁, λ₁)

Create a uniaxial hyperelastic data object from arrays of test data. The function returns a HyperelasticData object with the stresses and principal stretches. Currently, this assumes the material is incompressible.
"""
struct UniaxialHyperelasticData <: AbstractHyperelasticData
    s⃗
    λ⃗
    UniaxialHyperelasticData(s⃗, λ⃗) =
        let
            s⃗ = zip(s⃗)
            λ⃗ = zip(λ⃗, (λ⃗) .^ (-0.5), (λ⃗) .^ (-0.5))
            new(s⃗, λ⃗)
        end
end

struct BiaxialHyperelasticData <: AbstractHyperelasticData
    s⃗
    λ⃗
end

"""
BiaxialHyperelasticData(s₁, s₂, λ₁, λ₂)

Create a biaxial hyperelastic data object from arrays of test data. The function returns a HyperelasticData object with the stresses and principal stretches. Currently, this assumes the material is incompressible.
"""
function BiaxialHyperelasticData(s₁, s₂, λ₁, λ₂)
    s⃗ = zip(s₁, s₂)
    λ⃗ = zip(λ₁, λ₂, (λ₁ .* λ₂) .^ (-1))
    return BiaxialHyperelasticData(s⃗, λ⃗)
end


"""
---
HyperelasticProblem(data::HyperelasticData, model, u₀, ps; loss=L2DistLoss(), agg=AggMode.Mean(), cons=(x, p) -> [true], kwargs...)
---
Returns an `OptimizationProblem` for solving with GalacticOptim.jl. `data` is the hyperelastic experimental data, `model` is the strain energy density as a function of the parameters (i.e. `f(p) = W(p)(λ⃗)`). `ps` is any hyperparameters for the model (currently not supported). `loss` defines the loss function to be used in the optimization. Currently defaults to the ``L^2``-norm between the predicted and experimental data. `agg` defines the aggregration mode of the errors, defaults to the mean of the errors. `cons` define any constrain equations involving the parameters of the model. `kwargs` are passed to `OptimizationProblem`. To set parameter bounds, use the keywords `lb` and `ub` respectively.

"""
function HyperelasticProblem(data::AbstractHyperelasticData, ψ, u₀, ps; loss=L2DistLoss(), agg=AggMode.Mean(), cons=(x, p) -> [true], kwargs...)
    s = hcat(collect.(data.s⃗)...)

    stresses_provided = size(s, 1)

    function ŝ(p)
        s⃗ = NominalStressFunction(ψ, p)
        s₁₂₃ = @. s⃗(collect(data.λ⃗))
        return hcat(s₁₂₃...)[1:stresses_provided, :]
    end

    f(p, _) = [value(loss, s, ŝ(p), agg)]
    func = OptimizationFunction(f, Optimization.AutoForwardDiff(), cons=cons)
    return OptimizationProblem(func, u₀, ps; kwargs...)
end

"""
---
HyperelasticProblem(data::Vector{HyperelasticData}, model, u₀, ps; loss=L2DistLoss(), agg=AggMode.Mean(), cons=(x, p) -> [true], kwargs...)
---
Returns an `OptimizationProblem` for solving with GalacticOptim.jl. `data` is a vector of hyperelastic data and fits to all sets equally, `model` is the strain energy density as a function of the parameters (i.e. `f(p) = W(p)(λ⃗)`). `ps` is any hyperparameters for the model (currently not supported). `loss` defines the loss function to be used in the optimization. Currently defaults to the ``L^2``-norm between the predicted and experimental data. `agg` defines the aggregration mode of the errors, defaults to the mean of the errors. `cons` define any constrain equations involving the parameters of the model. `kwargs` are passed to `OptimizationProblem`. To set parameter bounds, use the keywords `lb` and `ub` respectively.

"""
function HyperelasticProblem(data::Vector{AbstractHyperelasticData}, model, u₀, ps; loss=L2DistLoss(), agg=AggMode.Mean(), cons=(x, p) -> [true], kwargs...)

    s = getindex.(vcat(collect.(zip(getfield.(data, :s⃗)...))...), 1) |> transpose
    λ = collect.(vcat(collect.(zip(getfield.(data, :λ⃗)...))...))

    stresses_provided = size(s, 1)

    function ŝ(p)
        s⃗ = NominalStressFunction(ψ, p)
        s₁₂₃ =  @. s⃗(collect(data.λ⃗))
        return hcat(s₁₂₃...)[1:stresses_provided, :]
    end

    f(p, _) = [value(loss, s, ŝ(p), agg)]
    func = OptimizationFunction(f, Optimization.AutoForwardDiff(), cons=cons)
    return OptimizationProblem(func, u₀, ps; kwargs...)
end

end
