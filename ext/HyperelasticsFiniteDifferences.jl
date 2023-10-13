module HyperelasticsFiniteDifferences

using FiniteDifferences
using Hyperelastics
using ADTypes

function Hyperelastics.∂ψ(
    ψ::Hyperelastics.AbstractHyperelasticModel{R},
    λ⃗::Vector{S},
    p,
    ad_type::ADTypes.AutoFiniteDifferences{T};
    kwargs...
) where {R, T,S}
    W(λ⃗) = StrainEnergyDensity(ψ, λ⃗, p)
    return grad(ad_type.fdm, W, λ⃗)[1]
end

end
