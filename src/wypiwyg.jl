export SussmanBathe
using DataInterpolations
using Quadrature
using QuadGK
using FiniteDifferences

function SussmanBathe(s⃗, λ⃗, k)
    e = log.(getindex.(λ⃗, 1))
    τ = DataInterpolations.CubicSpline(getindex.(s⃗, 1), e)
    w′(e) = 
        sum(
        i->τ((0.25^i)*e)+τ((-0.5)*((0.25)^i)*e), 
        0:k
        )
    probs = QuadratureProblem.(
        (u,p) -> w′(log(u)),
        [1.0, 1.0, 1.0],
        [2.0, 1/sqrt(2), 1/sqrt(2)],
    )
    map(prob -> solve(prob, QuadGKJL()), probs)
    # (λ⃗) -> sum(solve(remake(prob, lb = ones(3), ub = λ⃗), HCubatureJL(), reltol = 1e-8, abstol = 1e-8))
    (λ⃗) -> sum(
            map(i -> 
                solve(
                    remake(
                        probs[i], 
                        lb = 1.0,
                        ub = λ⃗[i]
                    ), 
                    QuadGKJL(),
                ).u,
                eachindex(probs)
            )
        )
end