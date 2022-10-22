

"""
First stretch invariant - Currently requires the addition of 5 times the machine precision to allow AD to work correctly

``I_1(\\vec{\\lambda}) = \\lambda_1^2+\\lambda_2^2+\\lambda_3^2 + 5\\varepsilon``
"""
ContinuumModels.I₁(λ⃗::AbstractVector) = sum(λ⃗ .^ 2)

"""
Second Stretch invariant

``I_2(\\vec{\\lambda}) = \\lambda_1^{-2}+\\lambda_2^{-2}+\\lambda_3^{-2}``
"""
ContinuumModels.I₂(λ⃗::AbstractVector) = sum(λ⃗ .^ (-2))

"""
Third Stretch invariant

``I_3(\\vec{\\lambda}) = (\\lambda_1\\lambda_2\\lambda_3)^2``
"""
ContinuumModels.I₃(λ⃗::AbstractVector) = prod(λ⃗)^2

"""
Volumetric Stretch

``J(\\vec{\\lambda}) = \\lambda_1\\lambda_2\\lambda_3``
"""
ContinuumModels.J(λ⃗::AbstractVector) = prod(λ⃗)
