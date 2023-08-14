module HyperelasticsEnzymeExt

using Enzyme
using Hyperelastics
using LinearAlgebra
using ADTypes
function Hyperelastics.∂ψ(
    ψ::Hyperelastics.AbstractHyperelasticModel,
    λ⃗::Vector{T},
    p,
    ad_type::AutoEnzyme;
    kwargs...,
) where {T}
    W(λ⃗) = StrainEnergyDensity(ψ, λ⃗, p)
    ∂W∂λ = collect(gradient(ad_type.mode, W, λ⃗, kwargs...))
    return ∂W∂λ
end

end
