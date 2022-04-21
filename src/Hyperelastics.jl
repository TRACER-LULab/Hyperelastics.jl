module Hyperelastics
using Reexport
using LossFunctions
using GalacticOptim
using ForwardDiff
using AbstractDifferentiation
using Tullio
using OMEinsum
# using Enzyme
# using Zygote
using ComponentArrays
# Write your package code here.
export HyperelasticData, uniaxial_data, biaxial_data, HyperelasticProblem
export s⃗̂, I₁, I₂, I₃, J
struct HyperelasticData
    s⃗
    λ⃗
end

function biaxial_data(s₁, s₂, λ₁, λ₂)
    s⃗ = zip(s₁, s₂)
    λ⃗ = zip(λ₁, λ₂, (λ₁ .* λ₂) .^ (-1))
    return HyperelasticData(s⃗, λ⃗)
end

function uniaxial_data(s₁, λ₁)
    s⃗ = zip(s₁)
    λ⃗ = zip(λ₁, (λ₁) .^ (-0.5), (λ₁) .^ (-0.5))
    return HyperelasticData(s⃗, λ⃗)
end

I₁(λ⃗) = sum(λ⃗ .^ 2) + 5eps(Float64)
I₂(λ⃗) = sum(λ⃗ .^ (-2)) + 5eps(Float64)
I₃(λ⃗) = prod(λ⃗)^2
J(λ⃗) = prod(λ⃗)

function s⃗̂(model, p, λ⃗; adb = AD.ForwardDiffBackend())
    W = model(p)
    σ₁₂₃ = map(x⃗ -> AD.gradient(adb, W, x⃗)[1].*x⃗,  λ⃗)
    σ̄₁₂₃ = map(x -> [x[1]-x[3], x[2]-x[3], x[3]-x[3]], σ₁₂₃)
    s₁₂₃ = map(x -> x[1]./ x[2], zip(σ̄₁₂₃, λ⃗))
end

include("hyperelastic_models.jl")
include("wypiwyg.jl")
function HyperelasticProblem(data::HyperelasticData, model, u₀, ps; loss=L2DistLoss(), agg=AggMode.Mean(), kwargs...)
    s = hcat(collect.(data.s⃗)...)
    stresses_provided = size(s, 1)
    s⃗(p) = s⃗̂(model, p, collect.(data.λ⃗))

    function ŝ(p)
        s₁₂₃ = s⃗(p)
        return hcat(s₁₂₃...)[1:stresses_provided, :]
    end

    f(p, _) = [value(loss, s, ŝ(p), agg)]
    func = OptimizationFunction(f, GalacticOptim.AutoForwardDiff())
    return OptimizationProblem(func, u₀, ps; kwargs...)
end
end
