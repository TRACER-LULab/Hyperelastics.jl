export SussmanBathe
using DataInterpolations
using QuadGK
using LossFunctions
"""
Sussman Bathe[^1]

Parameters: s⃗, λ⃗, k

Model: See Sussman-Bathe paper for description. For implementation, use s⃗̂(..., adb = AD.FiniteDifferencesBackend()).
See [link](https://github.com/JuliaDiff/ForwardDiff.jl/issues/510) and [link](https://github.com/SciML/Quadrature.jl/issues/56) for why Auto-diff does not work yet. Number of gausss points for SEF integration is specified by `k`

[^1]: > Sussman T, Bathe KJ. A model of incompressible isotropic hyperelastic material behavior using spline interpolations of tension–compression test data. Communications in numerical methods in engineering. 2009 Jan;25(1):53-63.
"""
function SussmanBathe(p)
    (; s⃗, λ⃗, k) = p
    s = DataInterpolations.CubicSpline(getindex.(s⃗, 1), getindex.(λ⃗, 1))
    function w′(λ)
        f(i) = s(λ^((1 / 4)^(i))) + s(λ^((-1 / 2) * ((1 / 4)^(i))))
        sum(f, 0:k)
    end
    w(x⃗) = sum(x -> quadgk(w′, 1.0, x)[1], x⃗)
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
function LatoreeMontans((;uniaxial, shear))
    if uniaxial.type == :T1
        w₁, w₃, w₁′, w₃′ = LatoreeMontansT1(uniaxial)
    elseif uniaxial.type == :T2A
        w₁, w₃, w₁′, w₃′ = LatoreeMontansT2A(uniaxial)
    elseif uniaxial.type == :T2B
        w₁, w₃, w₁′, w₃′ = LatoreeMontansT2A(uniaxial)
    end

    if shear.type == :PureShear
        w₁₃ = LatoreeMontansPS(shear)
    elseif shear.type == :SimpleShear
        w₁₃ = LatoreeMontansSS(shearm, w₁′, w₃′)
    end

    W(E₁, E₂, E₃, E₁₃, E₃₁, E₂₃, E₃₂) = w₁(E₁)+w₁(E₂)+w₃(E₃)+2w₁₃((E₁₃+E₃₁)/2)+2w₁₃((E₂₃_+E₃₂)/2)
end

function LatoreeMontansT1((; σ̃₁, Ẽ₁, Ẽ₂, k))
    # Integration Helper
    ∫(f, a, b) = quadgk(f, a, b)[1]
    # Intepolated Relationships from Stress-Strain Data
    σ₁ = DataInterpolations.CubicSpline(σ̃₁, Ẽ₁)
    E₂ = DataInterpolations.CubicSpline(Ẽ₂, Ẽ₁)
    E₃ = DataInterpolations.CubicSpline(-Ẽ₂ - Ẽ₁, Ẽ₁)
    # Find w₁
    # Equation 47
    function f₁(E₁, k)
        if k - 1 > 0
            f₁(E₂(E₁), k - 1)
        else
            return E₂(E₁)
        end
    end
    w₁′(E₁) = sum(σ₁ ∘ Base.Fix1(f₁, E₁), 0:k)
    w₁(E₁) = ∫(w₁′, 0, E₁)
    # Find w₃
    # From Equation 46
    W̃₃′ = w₁′.(E₂.(Ẽ₁))
    W₃′ = DataInterpolations.CubicSpline(w̃₃′, E₃.(Ẽ₁))
    w₃(E₃) = ∫(w₃′, 0, E₃)
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
        if k - 1 >= 0
            f₃(k - 1, E(E₃))
        else
            return E(E₃)
        end
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
            if k - 1 >= 0
                f₃(k - 1, E(E₃))
            else
                return E(E₃)
            end
        end
        w₃′(E₃) = sum(σ₃ ∘ Base.Fix2(f₃, E₃), 0:k)
        return value(L2DistLoss(), w₁′(E₂(Ẽ₁)), w₃(E₃(Ẽ₁)), AggMode.Sum())
    end
    optfunc = OptimizationFunction(f, Optimization.AutoForwardDiff())
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
        if k - 1 >= 0
            f₃(k - 1, E(E₃))
        else
            return E(E₃)
        end
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
function LatoreeMontansSS((;σ̃₁₃, γ̃), w₁′, w₃′)
    # σ₁₃ = DataInterpolations.CubicSpline(σ̃₁₃, γ̃)
    ψ̃ = @. 1/2*atan(2/γ̃)
    Ẽ₁₃ = @. -log(tan(ψ̃))*sin(2ψ̃)
    σ₁₃ = DataInterpolations.CubicSpline([σ̃₁₃; -σ̃₁₃], [Ẽ₁₃; Ẽ₁₃])
    ψ = DataInterpolations.CubicSpline([π/2-ψ̃, ψ̃], [-Ẽ₁₃, Ẽ₁₃])
    E₁(ψ) = -log(tan(ψ))*cos(2ψ)
    E₃(ψ) = log(tan(ψ))*cos(2ψ)
    w₁₃′(E₁₃) = (σ₁₃(E₁₃)-1/2*(w₁′(E₁(ψ(E₁₃))) - w₃′(E₃(ψ(E₁₃))))*(cos(2ψ(E₁₃))+E₁₃*sin(2ψ(E₁₃)))*sin(2ψ(E₁₃)))/(sin(2ψ(E₁₃))^2*(1-E₁(ψ(E₁₃))))
    w₁₃(E₁₃) = ∫(w₁₃′, 0, E₁₃)
end
