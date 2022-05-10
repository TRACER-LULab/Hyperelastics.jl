
export I₁, I₂, I₃, J, s⃗̂
"""
First stretch invariant - Currently requires the addition of 5 times the machine precision to allow AD to work correctly

``I_1 = \\lambda_1^2+\\lambda_2^2+\\lambda_3^2 + 5\\varepsilon``
"""
function I₁(λ⃗)
    sum(λ⃗ .^ 2) + 5eps(Float64)
end

"""
Second Stretch invariant

``I_2 = \\lambda_1^{-2}+\\lambda_2^{-2}+\\lambda_3^{-2}``
"""
function I₂(λ⃗)
    sum(λ⃗ .^ (-2)) + 5eps(Float64)
end

"""
Third Stretch invariant

``I_3 = (\\lambda_1\\lambda_\\lamdba_3)^2``
"""
function I₃(λ⃗)
    prod(λ⃗)^2
end

"""
Volumetric Stretch

``J = \\lambda_1\\lambda_2\\lambda_3``
"""
function J(λ⃗)
    prod(λ⃗)
end

function s⃗̂(model, p, λ⃗; adb=AD.ForwardDiffBackend())
    W = model(p)
    σ₁₂₃ = map(x⃗ -> AD.gradient(adb, W, x⃗)[1] .* x⃗, λ⃗)
    σ̄₁₂₃ = map(x -> [x[1] - x[3], x[2] - x[3], x[3] - x[3]], σ₁₂₃)
    s₁₂₃ = map(x -> x[1] ./ x[2], zip(σ̄₁₂₃, λ⃗))
    return s₁₂₃
end