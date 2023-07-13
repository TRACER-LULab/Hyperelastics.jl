module HyperelasticsForwardDiffExt

import ForwardDiff: gradient
using Hyperelastics
using ADTypes

function Hyperelastics.∂ψ(
    ψ::Hyperelastics.AbstractHyperelasticModel{R},
    λ⃗::Vector{T},
    p,
    ad_type::AutoForwardDiff;
    kwargs...,
) where {R,T}
    W(λ⃗) = StrainEnergyDensity(ψ, λ⃗, p)
    ∂W∂λ = gradient(W, λ⃗; kwargs...)
    return ∂W∂λ
end

end
