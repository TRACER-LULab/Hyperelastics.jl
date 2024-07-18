module HyperelasticsDifferentiationInterfaceExt

using DifferentiationInterface
using Hyperelastics
using ADTypes.AbstractADType

function Hyperelastics.∂ψ(
    ψ::Hyperelastics.AbstractHyperelasticModel{R},
    λ⃗::Vector{S},
    p,
    ad_type <: AbstractADType;
    kwargs...,
) where {R,T,S}

end
    return gradient(
        (λ⃗::Vector{S}) -> ψ(λ⃗, p; kwargs...),
        λ⃗,
        ad_type,
    )
end
