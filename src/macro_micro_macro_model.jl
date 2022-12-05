export MacroMicroMacro

struct MacroMicroMacro{T} <: AbstractDataDrivenHyperelasticModel
    Pchain::T
    r⃗::Vector{Vector{Float64}}
    w::Vector{Float64}
end

"""
`MacroMicroMacro(tests::Vector{AbstractHyperelasticTest}, Pchain, p₀;
                optimizer=LBFGS(),
                loss=L2DistLoss(),
                agg=AggMode.Mean())`

Model:
- See paper for description of data-driven model

Parameters:
- None

Fields:
- `tests`: Vector of Hyperelastic tests
- `Pchain(Pch, λch)`: Returns the interpolating function Pch(λch)
- `p₀``: Initial guess for discrete values of Pchain
- `optimizer`: Optimizer supported by [`Optimization.jl`](https://docs.sciml.ai/Optimization/stable/)
- `loss`: Loss on the differences between the predicted and experimental stresses. Uses a form from ['LossFunctions.jl`](https://juliaml.github.io/LossFunctions.jl/stable/)

> Amores VJ, Benítez JM, Montáns FJ. Data-driven, structure-based hyperelastic manifolds: A macro-micro-macro approach. arXiv preprint arXiv:1903.11545. 2019 Mar 27.
"""
function MacroMicroMacro(tests::Vector{AbstractHyperelasticTest}, Pchain, p₀; optimizer=LBFGS(), loss=L2DistLoss())
    λchain(λ⃗, r⃗) = dot(r⃗ .* r⃗, λ⃗)   # Eq. (1)
    ∂λchain∂λ(r⃗) = r⃗ .* r⃗         # In text below Eq. (8)
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
    λ⃗ = vcat(map(test -> test.data.λ, tests)...)
    @tullio λch[i, j] := λchain(λ⃗[i], r⃗[j])

    # Pre-compute Coefficients for quadrature
    @tullio wrᵢ²[k, i] := w[k] .* r⃗[k][i] .^ 2

    function calc_stresses(test, Pchain)
        λ⃗ = test.data.λ
        @tullio ∂W∂λ[i, j] := wrᵢ²[k, j] .* Pchain(λchain(λ⃗[i], r⃗[k])) # state i, stress j
        @tullio ŝ[i, j] := ∂W∂λ[j, i] - ∂W∂λ[j, 3] * λ⃗[j][3] / λ⃗[j][i] # state i, stress j
        return ŝ
    end

    function MSE(p, params)
        # Unpack Parameters
        (; s, tests, Pchain) = params
        # Create the Parameterized Pchain function
        _Pchain = Pchain(p)
        # Calculate the predicted stresses for each test
        ŝ = map(Base.Fix2(calc_stresses, _Pchain), tests)
        # Calculate the loss
        res = mean(map(i -> mean(abs, value(loss, s[i], ŝ[i][1:size(s[i], 1), :], AggMode.Mean(), ObsDim.First())), eachindex(s)))
        res += mean(p .< 0.0)
        return [res]
    end

    # Extract the stresses from the tesets
    s = map(x -> getfield.(x, :s), getfield.(tests, :data))
    s = map(Base.splat(hcat), s)

    # Solve the problem to find the optimal coefficients
    func = OptimizationFunction(
        MSE,
        Optimization.AutoForwardDiff()
    )
    prob = OptimizationProblem(
        func,
        p₀,
        (
            s=s,
            tests=tests,
            Pchain=Pchain,
        )
    )
    sol = solve(prob, optimizer)
    @show sol.retcode
    @assert sol.retcode == Optimization.ReturnCode.Success "Optimization of Parameters Failed"
    # Create the optimal PChain function
    pchain = Pchain(sol.u)
    MacroMicroMacro{typeof(pchain)}(pchain, r⃗, w)
end

function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::MacroMicroMacro, λ⃗::AbstractVector, p; kwargs...)
    f(x) = x[1] .* ψ.Pchain(sum(λ⃗ .* x[2] .^ 2)) .* x[2] .^ 2
    s = sum(map(f, zip(ψ.w, ψ.r⃗)))
    return s
end

function NonlinearContinua.CauchyStressTensor(ψ::MacroMicroMacro, λ⃗::AbstractVector, p; kwargs...)
    f(x) = x[1] .* ψ.Pchain(sum(λ⃗ .* x[2] .^ 2)) .* x[2] .^ 2
    s = sum(map(f, zip(ψ.w, ψ.r⃗)))
    σ = @. s * λ⃗
    return σ
end

function parameters(ψ::MacroMicroMacro)
    return nothing
end
