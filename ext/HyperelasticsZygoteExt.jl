module HyperelasticsZygoteExt

import Zygote: gradient
using Hyperelastics
using LinearAlgebra
using ADTypes

function Hyperelastics.∂ψ(ψ::Hyperelastics.AbstractHyperelasticModel, λ⃗::Vector{R}, p, ad_type::AutoZygote; kwargs...) where {R}
    W(λ⃗) = StrainEnergyDensity(ψ, λ⃗, p)
    ∂W∂λ = gradient(W, λ⃗, kwargs...)[1]
    return ∂W∂λ
end

# function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, F::AbstractMatrix, p, ad_type::AutoZygote; kwargs...)
#     σ = CauchyStressTensor(ψ, F, p, ad_type; kwargs...)
#     S = sqrt(det(F' * F)) * inv(F) * σ
#     return S
# end

# function NonlinearContinua.CauchyStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, λ⃗::AbstractVector, p, ad_type::AutoZygote; kwargs...)
#     S = SecondPiolaKirchoffStressTensor(ψ, λ⃗, p, ad_type; kwargs...)
#     σ = S .* λ⃗ ./ J(λ⃗)
#     return σ
# end

# function NonlinearContinua.CauchyStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, F::AbstractMatrix, p, ad_type::AutoZygote; kwargs...)
#     B = F * F'
#     a = eigvecs(B)'
#     B_prin = Diagonal(a * B * a')
#     V_prin = sqrt.(B_prin)
#     V = a' * V_prin * a
#     R = inv(V) * F
#     λ⃗ = sqrt.(diag(B_prin))
#     σ̂ = CauchyStressTensor(ψ, λ⃗, p, ad_type; kwargs...) |> Diagonal
#     σ = R * a' * σ̂ * a * R'
#     return σ
# end

end
