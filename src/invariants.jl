

"""
First stretch invariant - Currently requires the addition of 5 times the machine precision to allow AD to work correctly

``I_1(\\vec{\\lambda}) = \\lambda_1^2+\\lambda_2^2+\\lambda_3^2 + 5\\varepsilon``
"""
function I₁(λ⃗)
    sum(λ⃗ .^ 2)# + 10eps(Float64)
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
