module HyperelasticsEnzymeExt

# import NonlinearContinua
using Enzyme
using Hyperelastics
using LinearAlgebra

function Hyperelastics.∂ψ(
    ψ::Hyperelastics.AbstractHyperelasticModel,
    λ⃗::Vector{T},
    p,
    ad_type::AutoEnzyme;
    mode = Forward,
    kwargs...,
) where {T}
    W(λ⃗) = StrainEnergyDensity(ψ, λ⃗, p)
    ∂W∂λ = collect(gradient(mode, W, λ⃗, kwargs...))
    return ∂W∂λ
end

# function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, F::Matrix{T}, p, ad_type::AutoEnzyme; kwargs...) where {T}
#     σ = CauchyStressTensor(ψ, F, p, ad_type; kwargs...)
#     S = sqrt(det(F' * F)) * inv(F) * σ
#     return S
# end

# function NonlinearContinua.CauchyStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, λ⃗::Vector{T}, p, ad_type::AutoEnzyme; kwargs...) where {T}
#     S = SecondPiolaKirchoffStressTensor(ψ, λ⃗, p, ad_type; kwargs...)
#     σ = S .* λ⃗ ./ J(λ⃗)
#     return σ
# end

# function NonlinearContinua.CauchyStressTensor(ψ::Hyperelastics.AbstractHyperelasticModel, F::Matrix{T}, p, ad_type::AutoEnzyme; kwargs...) where {T}
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
