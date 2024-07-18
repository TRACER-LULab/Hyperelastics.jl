module HyperelasticsDifferentiationInterfaceExt

using DifferentiationInterface
using Hyperelastics
using ADTypes

function Hyperelastics.∂ψ(
    ψ::Hyperelastics.AbstractHyperelasticModel{R},
    λ⃗::Vector{T},
    p,
    ad_type::ADTypes.AbstractADType;
    kwargs...,
) where {R,T}
    W(λ⃗) = StrainEnergyDensity(ψ, λ⃗, p)
    ∂W∂λ = gradient(W, ad_type, λ⃗; kwargs...)
    return ∂W∂λ
end

end
