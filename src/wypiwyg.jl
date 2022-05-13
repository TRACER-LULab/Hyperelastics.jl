export SussmanBathe
using DataInterpolations
using Quadrature
using QuadGK

function SussmanBathe(p)
    (; s⃗, λ⃗, k) = p
    τ = DataInterpolations.CubicSpline(getindex.(s⃗, 1), getindex.(λ⃗, 1))
    w′(λ) = 
        sum(
        i -> τ(λ^(0.25^i))+τ(λ^((-0.5)*0.25^i)),
        0:k
        )
    probs = QuadratureProblem.(
        (u,p) -> w′(u),
        [1.0, 1.0, 1.0],
        [2.0, 1/sqrt(2), 1/sqrt(2)],
    )
    map(prob -> solve(prob, QuadGKJL()), probs)
    # (λ⃗) -> sum(solve(remake(prob, lb = ones(3), ub = λ⃗), HCubatureJL(), reltol = 1e-8, abstol = 1e-8))
    (λ⃗) -> sum(
            map(x -> 
                solve(
                    remake(
                        x[1], 
                        lb = 1.0,
                        ub = x[2]
                    ), 
                    QuadGKJL(),
                ).u,
                zip(probs, λ⃗)
            )
        )
end