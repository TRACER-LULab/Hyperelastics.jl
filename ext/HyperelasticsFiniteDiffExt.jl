module HyperelasticsFiniteDiffExt

import FiniteDiff: finite_difference_gradient
using Hyperelastics
using LinearAlgebra
using ADTypes

function Hyperelastics._SecondPiolaKirchoffStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel{T}, λ⃗::Vector{S}, p, ad_type::AutoFiniteDiff; return_type=eltype(λ⃗), kwargs...) where {T<:PrincipalValueForm, S}
    W(λ⃗) = StrainEnergyDensity(ψ, λ⃗, p)
    ∂W∂λ = finite_difference_gradient(W, λ⃗, ad_type.fdtype, return_type; kwargs...)
    return ∂W∂λ
end

end
