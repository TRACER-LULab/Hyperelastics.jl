export SussmanBathe
using DataInterpolations
# using Quadrature
using QuadGK
"""
Sussman Bathe

Parameters: s⃗, λ⃗, k

Optional Parameter: n 

Model: See Sussman-Bathe paper for description. For implementation, use s⃗̂(..., adb = AD.FiniteDifferencesBackend()). 
See [link](https://github.com/JuliaDiff/ForwardDiff.jl/issues/510) and [link](https://github.com/SciML/Quadrature.jl/issues/56) for why Auto-diff does not work yet. Number of gausss points for SEF integration is specified by `n`
"""
function SussmanBathe(p)
    (; s⃗, λ⃗, k) = p
    s = DataInterpolations.CubicSpline(getindex.(s⃗, 1), getindex.(λ⃗, 1))
    function w′(λ)
        f(i) = s(λ^((1 / 4)^(i))) + s(λ^((-1 / 2) * ((1 / 4)^(i))))
        sum(f, 0:k) 
    end
    w(x⃗) = sum(x -> quadgk(w′, 1.0, x)[1], x⃗)
    # probs = QuadratureProblem.(
    #     (u,p) -> w′(u),
    #     [1.0, 1.0, 1.0],
    #     [2.0, 1/sqrt(2), 1/sqrt(2)],
    # )
    # map(prob -> solve(prob, QuadGKJL(), reltol=1e-6), probs)
    # (λ⃗) -> sum(
    #         map(x -> 
    #             solve(
    #                 remake(
    #                     x[1], 
    #                     lb = convert(typeof(x[2]),1.0),
    #                     ub = x[2]
    #                 ), 
    #                 QuadGKJL(),
    #             ).u,
    #             zip(probs, λ⃗)
    #         )
    #     )
end