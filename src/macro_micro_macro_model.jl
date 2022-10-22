export MacroMicroMacro
"""
Macro-Micro-Macro Model [^1]

Model: See paper for description of data-driven model

Parameters:
`Pchain = f(p)(λchain)`: A higher-order function which returns the sets the parameters for the pchain function

`p₀`: Initial guess for `Pchain` parameters

`weights`: Weights the errors

`coeff_loss`: Loss function relating to the parameters

`solver`: The optimizer used to find the parameters for `Pchain` from Optimization.jl

`data`: Hyperelastic data for a biaxial test
---
[^1]: > Amores VJ, Benítez JM, Montáns FJ. Data-driven, structure-based hyperelastic manifolds: A macro-micro-macro approach. arXiv preprint arXiv:1903.11545. 2019 Mar 27.
"""
struct MacroMicroMacro <: AbstractDataDrivenHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::MacroMicroMacro, (; Pchain, p₀, weights, coeff_loss, solver, data))
    λchain(λ⃗, r⃗) = dot(r⃗ .* r⃗, λ⃗)   # Eq. (1)
    ∂λchain∂λ(λ⃗, r⃗) = r⃗ .* r⃗         # In text below Eq. (8)
    # 1b. Spherical quadrature from  `A micro-macro approach to rubber-like materials—Part I: the non-affine micro-sphere model of rubber elasticity`
    a = √(2) / 2
    b = 0.836095596749
    c = 0.387907304067
    r⃗ = [
        [0, 0, 1],
        [0, 1, 0],
        [1, 0, 0],
        [0, a, a],
        [0, -a, a],
        [a, 0, a],
        [-a, 0, a],
        [a, a, 0],
        [-a, a, 0],
        [b, c, c],
        [-b, c, c],
        [b, -c, c],
        [-b, -c, c],
        [c, b, c],
        [-c, b, c],
        [c, -b, c],
        [-c, -b, c],
        [c, c, b],
        [-c, c, b],
        [c, -c, b],
        [-c, -c, b],
    ]
    w1 = 0.0265214244093
    w2 = 0.0199301476312
    w3 = 0.0250712367487
    w = 2 .* [fill(w1, 3); fill(w2, 6); fill(w3, 12)] # Multiply by two since integration is over the half-sphere

    # Pre-compute all values of λch
    λ⃗ = collect(data.λ⃗)
    @tullio λch[i, j] := λchain(λ⃗[i], r⃗[j])

    # Pre-compute Coefficients for quadrature
    @tullio wrᵢ²[k, i] := w[k] .* r⃗[k][i] .^ 2

    function MSE(p, params)
        # Unpack Parameters
        data, λch, wrᵢ² = params
        s₁ = getindex.(data.s⃗, 1)
        s₂ = getindex.(data.s⃗, 2)

        # Create the Parameterized Pchain function
        _Pchain = Pchain(p)

        # Calculate the partial derivatives
        @tullio ∂W∂λᵢ[j, i] := wrᵢ²[k, i] .* _Pchain(λch[j, k])

        s̃₁ = @. ∂W∂λᵢ[:, 1] - ∂W∂λᵢ[:, 3] * getindex.(λ⃗, 3) / getindex.(λ⃗, 1)
        s̃₂ = @. ∂W∂λᵢ[:, 2] - ∂W∂λᵢ[:, 3] * getindex.(λ⃗, 3) / getindex.(λ⃗, 2)

        # Quadratic Error
        LossS1 = value(L2DistLoss(), s₁, s̃₁, AggMode.WeightedSum(weights))
        LossS2 = value(L2DistLoss(), s₂, s̃₂, AggMode.WeightedSum(weights))

        # Coefficient Loss
        LossC = coeff_loss(p)

        # Total Loss Term
        return [LossS1 + LossS2 + LossC]
    end

    # Solve the problem to find the optimal coefficients
    func = OptimizationFunction(MSE, Optimization.AutoForwardDiff())
    prob = OptimizationProblem(func, p₀, [data, λch, wrᵢ²])
    sol = solve(prob, solver, g_tol=1e-12)

    # Create the optimal PChain function
    pchain = Pchain(sol.u)

    # Create the SEF from the PChain Function
    Wch(λch) = quadgk(pchain, 1.0, λch)[1]
    W(λ⃗) = sum(x -> x[1] * Wch(λchain(λ⃗, x[2])), zip(w, r⃗))

    return W
end

# function NominalStressFunction(ψ::MacroMicroMacro, p)
# function ContinuumModels.StrainEnergyDensity(ψ::MacroMicroMacro, (; Pchain, p₀, weights, coeff_loss, solver, data))
#     λchain(λ⃗, r⃗) = dot(r⃗ .* r⃗, λ⃗)   # Eq. (1)
#     ∂λchain∂λ(λ⃗, r⃗) = r⃗ .* r⃗         # In text below Eq. (8)
#     # 1b. Spherical quadrature from  `A micro-macro approach to rubber-like materials—Part I: the non-affine micro-sphere model of rubber elasticity`
#     a = √(2) / 2
#     b = 0.836095596749
#     c = 0.387907304067
#     r⃗ = [
#         [0, 0, 1],
#         [0, 1, 0],
#         [1, 0, 0],
#         [0, a, a],
#         [0, -a, a],
#         [a, 0, a],
#         [-a, 0, a],
#         [a, a, 0],
#         [-a, a, 0],
#         [b, c, c],
#         [-b, c, c],
#         [b, -c, c],
#         [-b, -c, c],
#         [c, b, c],
#         [-c, b, c],
#         [c, -b, c],
#         [-c, -b, c],
#         [c, c, b],
#         [-c, c, b],
#         [c, -c, b],
#         [-c, -c, b],
#     ]
#     w1 = 0.0265214244093
#     w2 = 0.0199301476312
#     w3 = 0.0250712367487
#     w = 2 .* [fill(w1, 3); fill(w2, 6); fill(w3, 12)] # Multiply by two since integration is over the half-sphere

#     # Pre-compute all values of λch
#     λ⃗ = collect(data.λ⃗)
#     @tullio λch[i, j] := λchain(λ⃗[i], r⃗[j])

#     # Pre-compute Coefficients for quadrature
#     @tullio wrᵢ²[k, i] := w[k] .* r⃗[k][i] .^ 2

#     function MSE(p, params)
#         # Unpack Parameters
#         data, λch, wrᵢ² = params
#         s₁ = getindex.(data.s⃗, 1)
#         s₂ = getindex.(data.s⃗, 2)

#         # Create the Parameterized Pchain function
#         _Pchain = Pchain(p)

#         # Calculate the partial derivatives
#         @tullio ∂W∂λᵢ[j, i] := wrᵢ²[k, i] .* _Pchain(λch[j, k])

#         s̃₁ = @. ∂W∂λᵢ[:, 1] - ∂W∂λᵢ[:, 3] * getindex.(λ⃗, 3) / getindex.(λ⃗, 1)
#         s̃₂ = @. ∂W∂λᵢ[:, 2] - ∂W∂λᵢ[:, 3] * getindex.(λ⃗, 3) / getindex.(λ⃗, 2)

#         # Quadratic Error
#         LossS1 = value(L2DistLoss(), s₁, s̃₁, AggMode.WeightedSum(weights))
#         LossS2 = value(L2DistLoss(), s₂, s̃₂, AggMode.WeightedSum(weights))

#         # Coefficient Loss
#         LossC = coeff_loss(p)

#         # Total Loss Term
#         return [LossS1 + LossS2 + LossC]
#     end

#     # Solve the problem to find the optimal coefficients
#     func = OptimizationFunction(MSE, Optimization.AutoForwardDiff())
#     prob = OptimizationProblem(func, p₀, [data, λch, wrᵢ²])
#     sol = solve(prob, solver, g_tol=1e-12)

#     # Create the optimal PChain function
#     pchain = Pchain(sol.u)
#     function s⃗(λ⃗)
#         @tullio λch[i] := λchain(λ⃗, r⃗[i])
#         @tullio ∂W∂λᵢ[i] := wrᵢ²[j, i] .* Pchain(λch[j])
#         return @. ∂W∂λᵢ - ∂W∂λᵢ[3] * λ⃗[3] / λ⃗[1]
#     end
#     return s⃗
# end
# end

function parameters(ψ::MacroMicroMacro)
    return (:Pchain, :p₀, :weights, :coeff_loss, :solver, :data)
end
