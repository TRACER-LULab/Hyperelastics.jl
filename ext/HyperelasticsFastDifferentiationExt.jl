module HyperelasticsFastDifferentiationExt

using FastDifferentiation
import Hyperelastics
using LinearAlgebra
import ADTypes

function Hyperelastics.∂ψ(
    ψ::Hyperelastics.AbstractHyperelasticModel{T},
    λ⃗::Vector{S},
    p,
    ad_type;
    fast_diff = true,
    kwargs...,
) where {T,S}
    # @info "NOTE THIS ONLY RETURNS AN EXECUTABLE FUNCTION FOR GETTING THE VALUES REQUIRED"

    λ = make_variables(:λ, 3)
    W = Hyperelastics.StrainEnergyDensity(ψ, λ, p)
    ∂W∂λ = jacobian([W], λ)
    f = make_function(∂W∂λ, λ)

    return f(λ⃗) |> vec
end

end
