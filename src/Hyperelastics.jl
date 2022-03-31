module Hyperelastics
using Reexport
using Symbolics, ModelingToolkit, LossFunctions, GalacticOptim
using ForwardDiff
using LabelledArrays
# Write your package code here.
export HyperelasticData, uniaxial_data, biaxial_data, HyperelasticProblem
export NeoHookean, MooneyRivlin, Gent
struct HyperelasticData
    s⃗
    λ⃗
end

function biaxial_data(s₁, s₂, λ₁, λ₂)
    s⃗ = [s₁; s₂]>
    λ⃗ = hcat([λ₁, λ₂]...)
    return HyperelasticData(s⃗, λ⃗)
end

function uniaxial_data(s₁, λ₁)
    s⃗ = s₁
    λ⃗ = zip(λ₁, (λ₁) .^ (-0.5), (λ₁) .^ (-0.5))
    return HyperelasticData(s⃗, λ⃗)
end

I₁(λ⃗) = sum(λ⃗ .^ 2)
I₂(λ⃗) = sum(λ⃗ .^ -2)
include("hyperelastic_models.jl")

function HyperelasticProblem(data::HyperelasticData, model, u₀, ps; loss=L2DistLoss(), agg=AggMode.Mean(), kwargs...)
    function s(u)
        W = model(u)
        s⃗ = ForwardDiff.gradient.(W, collect.(data.λ⃗))
        λ₁ = getindex.(collect.(data.λ⃗), 1)
        λ₃ = getindex.(collect.(data.λ⃗), 3)
        s₁ = getindex.(collect.(s⃗), 1)
        s₃ = getindex.(collect.(s⃗), 3)
        σ11 = λ₁.*s₁
        σ33 = λ₃.*s₃
        (σ11.-σ33)./λ₁
    end
    f(u, p) = [value(loss, Vector(data.s⃗), s(u), agg)]
    func = OptimizationFunction(f, GalacticOptim.AutoForwardDiff())
    return OptimizationProblem(func, u₀, ps, kwargs...)
end

end
