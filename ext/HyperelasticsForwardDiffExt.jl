module HyperelasticsForwardDiffExt

import ForwardDiff: gradient

using Hyperelastics
using LinearAlgebra
using ADTypes

function Hyperelastics._SecondPiolaKirchoffStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel{R}, λ⃗::Vector{T}, p, ad_type::AutoForwardDiff; kwargs...) where {R<:PrincipalValueForm,T}
    W(λ⃗) = StrainEnergyDensity(ψ, λ⃗, p)
    ∂W∂λ = gradient(W, λ⃗; kwargs...)
    return ∂W∂λ
end

end
