module Hyperelastics
using Reexport
using Symbolics, ModelingToolkit, LossFunctions, GalacticOptim
using ForwardDiff
using LabelledArrays
# Write your package code here.
export HyperelasticData, uniaxial_data, biaxial_data, HyperelasticProblem
export NeoHookian, MooneyRivlin, Gent
struct HyperelasticData
    s⃗
    λ⃗
end

function biaxial_data(s₁, s₂, λ₁, λ₂)
    s⃗ = [s₁; s₂]
    λ⃗ = hcat([λ₁, λ₂]...)
    return HyperelasticData(s⃗, λ⃗)
end

function uniaxial_data(s₁, λ₁)
    s⃗ = s₁
    λ⃗ = zip(λ₁, (λ₁) .^ (-0.5), (λ₁) .^ (-0.5))
    return HyperelasticData(s⃗, λ⃗)
end

I₁(λ⃗) = sum(λ⃗.^2)
include("hyperelastic_models.jl")

function HyperelasticProblem(data::HyperelasticData, model, u₀, ps; loss=L2DistLoss(), agg=AggMode.Mean(), kwargs...)
    function s(u)
        W = model(u)
        s⃗ = ForwardDiff.gradient.(W, collect.(data.λ⃗))
        getindex.(s⃗, 1) - getindex.(s⃗, 3)
    end
    f(u, p) = [value(loss, Vector(data.s⃗), s(u), agg)]
    func = OptimizationFunction(f, GalacticOptim.AutoForwardDiff())
    return OptimizationProblem(func, u₀, ps, lb = LVector(μ = 0.0, Jₘ = 0.2178580576229847), ub = LVector(μ = Inf, Jₘ = Inf), kwargs...)
end

end
