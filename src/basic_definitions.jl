
"""
First stretch invariant - Currently requires the addition of 5 times the machine precision to allow AD to work correctly

``I_1(\\vec{\\lambda}) = \\lambda_1^2+\\lambda_2^2+\\lambda_3^2 + 5\\varepsilon``
"""
function I₁(λ⃗)
    sum(λ⃗ .^ 2)+10eps(Float64)
end

"""
Second Stretch invariant

``I_2(\\vec{\\lambda}) = \\lambda_1^{-2}+\\lambda_2^{-2}+\\lambda_3^{-2}``
"""
function I₂(λ⃗)
    sum(λ⃗ .^ (-2))
end

"""
Third Stretch invariant

``I_3(\\vec{\\lambda}) = (\\lambda_1\\lambda_2\\lambda_3)^2``
"""
function I₃(λ⃗)
    prod(λ⃗)^2
end

"""
Volumetric Stretch

``J(\\vec{\\lambda}) = \\lambda_1\\lambda_2\\lambda_3``
"""
function J(λ⃗)
    prod(λ⃗)
end

"""
s⃗̂(model, λ⃗; adb=AD.ForwardDiffBackend())

Return nominal stress predicted by a `model` at stretches, `λ⃗`, using differentiation mode, `adb` from AbstractDifferentiation.jl. Defaults to using ForwardDiff.jl for AD.
"""
# function s⃗̂(W, λ⃗; adb=AD.ForwardDiffBackend())
#     σ₁₂₃ = map(x⃗ -> AD.gradient(adb, W, x⃗)[1] .* x⃗, λ⃗)
#     σ̄₁₂₃ = map(x -> [x[1] - x[3], x[2] - x[3], x[3] - x[3]], σ₁₂₃)
#     s₁₂₃ = map(x -> x[1][1:3] ./ x[2][1:3], zip(σ̄₁₂₃, λ⃗))
#     return s₁₂₃
# end

function s⃗̂(W, λ⃗; adb = AD.ForwardDiffBackend())
    ∂W∂λᵢ = map(λ⃗ᵢ -> AD.gradient(adb, W, λ⃗ᵢ)[1], λ⃗)
    s⃗ᵢ = ∂W∂λᵢ .- getindex.(∂W∂λᵢ, 3) .* getindex.(λ⃗, 3) .* map(λ⃗ᵢ -> 1 ./ λ⃗ᵢ,  λ⃗)
end
