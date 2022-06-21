using DataInterpolations
using QuadGK
using LinearAlgebra
using LossFunctions
using NonlinearSolve
export SussmanBathe, LatoreeMontans, StabilitySmoothedSussmanBathe
"""
Sussman Bathe[^1]

Parameters: σ⃗, λ⃗, k

Model: See Sussman-Bathe paper for description. For implementation, use s⃗̂(..., adb = AD.FiniteDifferencesBackend()).
See [link](https://github.com/JuliaDiff/ForwardDiff.jl/issues/510) and [link](https://github.com/SciML/Quadrature.jl/issues/56) for why Auto-diff does not work yet. Number of gausss points for SEF integration is specified by `k`

[^1]: > Sussman T, Bathe KJ. A model of incompressible isotropic hyperelastic material behavior using spline interpolations of tension–compression test data. Communications in numerical methods in engineering. 2009 Jan;25(1):53-63.
"""
function SussmanBathe(p)
    (; s⃗, λ⃗, k) = p
    ∫(f, a, b) = quadgk(f, a, b, rtol=1e-8)[1]
    σ = DataInterpolations.CubicSpline(getindex.(s⃗, 1) .* getindex.(λ⃗, 1), getindex.(λ⃗, 1))
    f(i, λ) = (σ(λ^((4.0)^(-i))) + σ(λ^(-0.5 * (4.0^(-i))))) / λ
    function w′(λ)
        sum(Base.Fix2(f, λ), 0:k)
    end
    w(x⃗) = sum(x -> ∫(w′, 1.0, x), x⃗)
end

"""
Smoothed Sussman Bathe[^1]

Model: See paper for smoothing description.

Parameters: s⃗

[^1]: > Latorre M, Montáns FJ. Experimental data reduction for hyperelasticity. Computers & Structures. 2020 May 1;232:105919.
"""
function StabilitySmoothedSussmanBathe(p)
    (; s⃗, λ⃗, k, n, q) = p
    # Integration
    ∫(f, a, b) = quadgk(f, a, b)[1]

    #  Definitions
    Ê = log.(getindex.(λ⃗, 1))
    Emin, Emax = extrema(Ê)
    σ̂ = getindex.(s⃗, 1) .* getindex.(λ⃗, 1)
    m = n + k + 1
    # Create the Basis Functions
    t⃗ = range(extrema(Ê)..., length=m)
    N⃗(E) = map(i->N(E, k=k, i=i, t⃗=t⃗), 1:n)
    N̂(E⃗) = vcat(N⃗.(E⃗)...)

    # σ(E) Fitting
    σ(E, B) = N⃗(E) * B
    σ⃗(E⃗, B) = N̂(E⃗) * B
    σ′(E, B) = ForwardDiff.derivative(Base.Fix2(σ, B), E)
    σ′′(E, B) = ForwardDiff.derivative(Base.Fix2(σ′, B), E)
    σ′′2(E, B) = σ′′(E, B)^2

    # Regression Function
    f̂(E⃗, B; Wr=I) = 1 / 2 / length(E⃗) * (σ⃗(E⃗, B) - σ̂)' * Wr * (σ⃗(E⃗, B) - σ̂)

    # Curvature Based Smoothing
    f̃_curve(B) = (∫(Base.Fix2(σ′′2, B), Emin, Emax)) / (Emax - Emin)

    # Stability Based Smoothing
    Cp(E, B) = σ′(E, B) - σ(E, B)
    Cm(E, B) = 2σ′(E, B) - σ(E, B)

    dCpdE(E, B) = ForwardDiff.derivative(Base.Fix2(Cp, B), E)
    dSpdE(E, B) = exp(-2E) * (-2Cp(E, B) + dCpdE(E, B))

    dCmdE(E, B) = ForwardDiff.derivative(Base.Fix2(Cp, B), E)
    dSmdE(E, B) = exp(E) * (Cm(E, B) + dCmdE(E, B))

    f̃_stability(B) = (∫(x -> dSmdE(x, B)^2, Emin, 0) + ∫(x -> dSpdE(x, B)^2, 0, Emax)) / (Emax - Emin)

    # Smoothing Function
    f(B, _) = [(1-q)*f̂(Ê, B)]#+ q * f̃_curve(B)]
    # f(B, _) = [(1 - q) * f̂(σ⃗(Ê, B), B) + q * f̃_stability(B)]
    # Optimize the Control Point
    func = OptimizationFunction(f, Optimization.AutoForwardDiff())
    prob = OptimizationProblem(func, randn(n), [])
    sol = solve(prob, LBFGS())

    # Return the σ(E) for the optimal B
    return Base.Fix2(σ, sol.u)
end

function N(t; k, i, t⃗)
    if k ≥ 2
        return (t - t⃗[i]) / (t⃗[i+k-1] - t⃗[i]) * N(t, k=k - 1, i=i, t⃗=t⃗) + (t⃗[i+k] - t) / (t⃗[i+k] - t⃗[i+1]) * N(t, k=k - 1, i=i + 1, t⃗=t⃗)
    elseif t⃗[i] ≤ t < t⃗[i+1]
        return 1.0
    else
        return 0.0
    end
end
"""
Latoree Montans [^1] - Incompressible Transversely isotropic

Parameters:

Optional Parameters:

Model:
Type I - :TypeI - One Uniaxial Test in a transversely isotropic direction
- Parameters: σ̃₁, Ẽ₁, Ẽ₂, k
Type 2 - :TypeII - Two independent uniaxial tests in the isotropic and anisotropic directions where 3 is the preferred direction for a uniaxial tesnsion compression test.

[^1]: > Latorre M, Montáns FJ. Extension of the Sussman–Bathe spline-based hyperelastic model to incompressible transversely isotropic materials. Computers & Structures. 2013 Jun 1;122:13-26.
"""
function LatoreeMontans(uniaxial, shear)
    if uniaxial.type == :T1
        w₁, w₃, w₁′, w₃′ = LatoreeMontansT1(uniaxial)
    elseif uniaxial.type == :T2A
        w₁, w₃, w₁′, w₃′ = LatoreeMontansT2A(uniaxial)
    elseif uniaxial.type == :T2B
        w₁, w₃, w₁′, w₃′ = LatoreeMontansT2B(uniaxial)
    end
    if shear.type == :PureShear
        w₁₃ = LatoreeMontansPS(shear)
    elseif shear.type == :SimpleShear
        w₁₃ = LatoreeMontansSS(shear, w₁′, w₃′)
    else
        w₁₃(x) = 0.0
    end

    function W(E)
        w₁(E[1]) + w₁(E[2]) + w₁(E[3]) + 2w₁₃((E[4] + E[5]) / 2) + 2w₁₃((E[6] + E[7]) / 2)
    end
end

function LatoreeMontansT1((; σ̃₁, Ẽ₁, Ẽ₂, k))
    # Integration Helper
    ∫(f, a, b) = quadgk(f, a, b)[1]
    # Intepolated Relationships from Stress-Strain Data
    σ₁ = DataInterpolations.CubicSpline(σ̃₁, Ẽ₁)
    E₂ = DataInterpolations.CubicSpline(Ẽ₂, Ẽ₁)
    Ẽ₃ = 0.0 .- Ẽ₂ .- Ẽ₁
    E₃ = DataInterpolations.CubicSpline(Ẽ₃, Ẽ₁)
    # Find w₁
    # Equation 47
    function f₁(E₁, k)
        #k = 0
        res = E₁
        for _ in 1:k
            res = E₂(res)
        end
        return res
    end
    display(f₁(1.0, 4))
    w₁′(E₁) = sum(i -> σ₁(f₁(E₁, i)), 0:k)
    w₁(E₁) = ∫(w₁′, 0.0, E₁)
    # Find w₃
    # From Equation 46
    E₂̄ = E₂.(Ẽ₁)
    w̃₃′ = w₁′.(E₂̄)
    E₃̄ = E₃.(Ẽ₁)
    w₃′ = DataInterpolations.CubicSpline(w̃₃′, E₃̄)
    w₃(E₃) = ∫(w₃′, 0.0, E₃)
    return w₁, w₃, w₁′, w₃′
end
function LatoreeMontansT2A((; σ̃₃, Ẽ₁, Ẽ₂, Ẽ₃, k))
    # Integration
    ∫(f, a, b) = quadgk(f, a, b)[1]
    # interpolations
    σ₃ = DataInterpolations.CubicSpline(σ̃₃, Ẽ₃)
    E₂ = DataInterpolations.CubicSpline(Ẽ₂, Ẽ₁)
    E₃ = DataInterpolations.CubicSpline(-Ẽ₂ - Ẽ₁, Ẽ₂)
    # Find w₃′
    E(E₃) = E₃(-E₃ / 2)
    function f₃(k, E₃)
        #k = 0
        res = E₃
        for _ in 1:k
            res = E(E₃)
        end
        return res
    end
    w₃′(E₃) = sum(σ₃ ∘ Base.Fix2(f₃, E₃), 0:k)
    w₃(E₃) = ∫(w₃′, 0, E₃)
    # Find w₁′
    w₁′ = DataInterpolations.CubicSpline(w₃′.(E₃(E₂)), E₂)
    w₁(E₁) = ∫(w₁′, 0, E₁)
    return w₁, w₃, w₁′, w₃′
end
function LatoreeMontansT2B((; σ̃₁, σ̃₃, Ẽ₃, Ẽ₁, Ê₂, k, p₀, solver))
    # Ê₂(E₁, p) = ....
    # integration Helper
    ∫(f, a, b) = quadgk(f, a, b)[1]
    # Data interpolations
    σ₁ = DataInterpolations.CubicSpline(σ̃₁, Ẽ₁)
    σ₃ = DataInterpolations.CubicSpline(σ̃₃, Ẽ₃)
    # Function to optimize parameters of Ê₂
    function f(p, _)
        Ẽ₂ = map(Base.Fix2(Ê₂, p), Ẽ₁)
        E₂ = DataInterpolations.CubicSpline(Ẽ₂, Ẽ₁)
        function f₁(E₁, k)
            # k = 0
            res = E₁
            for _ in 1:k
                res = E₂(res)
            end
            return res
        end
        w₁′(E₁) = sum(i -> σ₁(f₁(E₁, i)), 0:k)
        E₃ = DataInterpolations.CubicSpline(-Ẽ₂ - Ẽ₁, Ẽ₂)
        # Find w₃′
        E(E₃) = E₃(-E₃ / 2)
        function f₃(k, E₃)
            res = E₃
            for _ in 1:k
                res = E(E₃)
            end
            return res
        end
        w₃′(E₃) = sum(i -> σ₃(f₃(i, E₃)), 0:k)
        return value(L2DistLoss(), w₁′(E₂(Ẽ₁)), w₃(E₃(Ẽ₁)), AggMode.Sum())
    end
    optfunc = OptimizationFunction(f, Optimization.AutoFiniteDiff())
    optprob = OptimizationProblem(optfunc, p₀, [])
    sol = solve(optprob, solver)
    E₂(E₁) = Base.Fix2(Ê₂, sol.u)
    Ẽ₂ = map(E₂, Ẽ₁)
    function f₁(E₁, k)
        if k - 1 > 0
            f₁(E₂(E₁), k - 1)
        else
            return E₂(E₁)
        end
    end
    w₁′(E₁) = sum(σ₁ ∘ Base.Fix1(f₁, E₁), 0:k)
    E₃ = DataInterpolations.CubicSpline(-Ẽ₂ - Ẽ₁, Ẽ₂)
    # Find w₃′
    E(E₃) = E₃(-E₃ / 2)
    function f₃(k, E₃)
        res = E(E₃)
        for _ in 1:k
            res = E(E₃)
        end
        return res
    end
    w₃′(E₃) = sum(σ₃ ∘ Base.Fix2(f₃, E₃), 0:k)
    w₁(E₁) = ∫(w₁′, 0, E₁)
    w₃(E₃) = ∫(w₃′, 0, E₃)
    return w₁, w₃, w₁′, w₃′
end
function LatoreeMontansPS((; σ̃, Ẽ))
    σ = DataInterpolations.CubicSpline([-σ̃; σ̃], [-Ẽ; Ẽ])
    w₁₃′(E₁₃) = σ(E₁₃)
    w₁₃(E₁₃) = ∫(w₁₃′, 0, E₁₃)
end
function LatoreeMontansSS((; σ̃₁₃, γ̃), w₁′, w₃′)
    # σ₁₃ = DataInterpolations.CubicSpline(σ̃₁₃, γ̃)
    ψ̃ = @. 1 / 2 * atan(2 / γ̃)
    Ẽ₁₃ = @. -log(tan(ψ̃)) * sin(2ψ̃)
    σ₁₃ = DataInterpolations.CubicSpline([σ̃₁₃; -σ̃₁₃], [Ẽ₁₃; Ẽ₁₃])
    ψ = DataInterpolations.CubicSpline([π / 2 - ψ̃, ψ̃], [-Ẽ₁₃, Ẽ₁₃])
    E₁(ψ) = -log(tan(ψ)) * cos(2ψ)
    E₃(ψ) = log(tan(ψ)) * cos(2ψ)
    w₁₃′(E₁₃) = (σ₁₃(E₁₃) - 1 / 2 * (w₁′(E₁(ψ(E₁₃))) - w₃′(E₃(ψ(E₁₃)))) * (cos(2ψ(E₁₃)) + E₁₃ * sin(2ψ(E₁₃))) * sin(2ψ(E₁₃))) / (sin(2ψ(E₁₃))^2 * (1 - E₁(ψ(E₁₃))))
    w₁₃(E₁₃) = ∫(w₁₃′, 0, E₁₃)
end
"""
Macro-Micro-Macro Model [^1]

Model: See paper for description of data-driven model

Parameters: s₁, λ₁ (uniaxial test data)+
`nv`: number of B-spline vertices
`nch`: number of chain in a sphere
`nq`: number of spherical quadrature points
---
[^1]: > Amores VJ, Benítez JM, Montáns FJ. Data-driven, structure-based hyperelastic manifolds: A macro-micro-macro approach. arXiv preprint arXiv:1903.11545. 2019 Mar 27.
"""
function MacroMicroMacro((; nv, nq))
    # 1. Determine Pch(λch) from uniaxial tensile test
    # 1a. λch => λ⃗, r⃗ Eq. 1
    λch(λ⃗, r⃗) = dot(λ⃗, r⃗ .^ 2)
    # 1b. Spherical quadrature from  `A micro-macro approach to rubber-like materials—Part I: the non-affine micro-sphere model of rubber elasticity`
    r = [
        [0, 0, 1],
        [0, 1, 0],
        [1, 0, 0],
        [0, 1 / √2, 1 / √2],
        [0, -1 / √2, 1 / √2],
        [1 / √2, 0, 1 / √2],
        [-1 / √2, 0, 1 / √2],
        [1 / √2, 1 / √2, 0],
        [-1 / √2, 1 / √2, 0],
        [0.836095596749, 0.387907304067, 0.387907304067],
        [-0.836095596749, 0.387907304067, 0.387907304067],
        [0.836095596749, -0.387907304067, 0.387907304067],
        [-0.836095596749, -0.387907304067, 0.387907304067],
        [0.387907304067, 0.836095596749, 0.387907304067],
        [-0.387907304067, 0.836095596749, 0.387907304067],
        [0.387907304067, -0.836095596749, 0.387907304067],
        [-0.387907304067, -0.836095596749, 0.387907304067],
        [0.387907304067, 0.387907304067, 0.836095596749],
        [-0.387907304067, 0.387907304067, 0.836095596749],
        [0.387907304067, -0.387907304067, 0.836095596749],
        [-0.387907304067, -0.387907304067, 0.836095596749],
    ]
    w = [
        0.0265214244093,
        0.0265214244093,
        0.0265214244093,
        0.0199301476312,
        0.0199301476312,
        0.0199301476312,
        0.0199301476312,
        0.0199301476312,
        0.0199301476312,
        0.0250712367487,
        0.0250712367487,
        0.0250712367487,
        0.0250712367487,
        0.0250712367487,
        0.0250712367487,
        0.0250712367487,
        0.0250712367487,
        0.0250712367487,
        0.0250712367487,
        0.0250712367487,
        0.0250712367487
    ]
    # 1c. BSpline Function
    function N(k, i, t, t⃗)
        if k ≥ 2
            (t - t⃗[i]) / (t[i+k-1] - t[i]) * N(k - 1, i, t, t⃗) + (t[i+k] - t) / (t[i+k] - t[i+1]) * N(k - 1, i + 2, t, t⃗)
        elseif t[i] ≤ t < t[i]
            1
        else
            0
        end
    end
    # 1d.
end

function SphericalQuadrature(nq)
    function f(x, p)
        # x = reshape(x, 4, :)
        r1 = sum(i -> x[i:i+2] * x[i+3], 1:nq)
        r2 = sum(i -> (x[i:i+2] * (x[i:i+2]')) * x[i+3], 1:nq)
        return [reshape(r1, :, 1); reshape(r2, :, 1)]
        # r1 = sum(i->x.r[i]*x.w[i], 1:nq)
        # r2 = sum(i->(x.r[i]*x.r[i]')*w[i], 1:nq)-1/3*I(3)
        # return [reshape(r2, :, 1); reshape(r1,:, 1)]
    end
    # r = map(x->ComponentArray(
    #         x = randn(),
    #         y = randn(),
    #         z = randn()
    #         ), 1:nq)
    # w = rand(nq)
    # x₀ = ComponentArray(r=r, w=w)
    # display(f(x₀, []))
    x₀ = randn(nq * 4)
    display(f(x₀, []))
    func = NonlinearFunction(f)
    np = NonlinearProblem(func, x₀, [])
    sol = solve(np, NewtonRaphson())
end
