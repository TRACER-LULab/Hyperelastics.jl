"""
---
HyperelasticProblem(data::HyperelasticData, model, u₀, ps; loss=L2DistLoss(), agg=AggMode.Mean(), cons=(x, p) -> [true], kwargs...)
---
Returns an `OptimizationProblem` for solving with GalacticOptim.jl. `data` is the hyperelastic experimental data, `model` is the strain energy density as a function of the parameters (i.e. `f(p) = W(p)(λ⃗)`). `ps` is any hyperparameters for the model (currently not supported). `loss` defines the loss function to be used in the optimization. Currently defaults to the ``L^2``-norm between the predicted and experimental data. `agg` defines the aggregration mode of the errors, defaults to the mean of the errors. `cons` define any constrain equations involving the parameters of the model. `kwargs` are passed to `OptimizationProblem`. To set parameter bounds, use the keywords `lb` and `ub` respectively.

"""
function HyperelasticProblem(
    data::AbstractHyperelasticData,
    ψ::AbstractHyperelasticModel,
    u₀,
    ps;
    loss=L2DistLoss(),
    agg=AggMode.Mean(),
    ad=Optimization.AutoForwardDiff(),
    kwargs...
)

    stress_provided = length(data.s⃗[1])

    function ŝ(p)
        s⃗ = map(x -> SecondPiolaKirchoffStressTensor(ψ, x, p), data.λ⃗)
        λ⃗ = data.λ⃗

        @tullio σ⃗[i, j] := s⃗[i][j] .* λ⃗[i][j]
        @tullio Δσ[i, j, k] := σ⃗[i, j] - σ⃗[i, k]
        @tullio Δs[i, j, k] := Δσ[i, j, k] / λ⃗[i][j]

        Δs₁₃ = Δs[:, 1, 3]
        Δs₂₃ = Δs[:, 2, 3]
        Δs₁₂ = Δs[:, 1, 2]

        res = map(x -> [x[1], x[2], x[3]], zip(Δs₁₃, Δs₂₃, Δs₁₂))
        res = hcat(res...)
        return res[1:stress_provided, :]
    end

    f(p, _) = [value(loss, hcat(data.s⃗...), ŝ(p), agg)]

    cons = constraints(ψ, data)
    lb, ub = parameter_bounds(ψ, data)
    func = OptimizationFunction(f, ad)
    prob = OptimizationProblem(func, u₀, ps)
    # Check for Constraints
    if !isnothing(cons)
        # println("Has Constraints")
        num_cons = length(cons(u₀, ps))
        func = OptimizationFunction(f, ad, cons=cons)
        prob = OptimizationProblem(func, u₀, ps, lcons=zeros(num_cons))
    end
    # Check for Bounds
    if !isnothing(lb) || !isnothing(ub)
        # println("Has Bounds")
        ax = Axis(Hyperelastics.parameters(ψ))
        if !isnothing(lb) && !isnothing(ub)
            lb = LVector(lb)

        elseif !isnothing(lb)
            lb = ComponentVector(lb)
            ub = ComponentVector(ones(length(lb)) * Inf, ax)
        elseif !isnothing(ub)
            ub = ComponentVector(ub)
            lb = ComponentVector(ones(length(ub)) * -Inf, ax)
        end
        if isnothing(prob.lcons)
            prob = OptimizationProblem(func, u₀, ps, lb=lb, ub=ub)
        else
            num_cons = length(cons(u₀, ps))
            prob = OptimizationProblem(func, u₀, ps, lcons=zeros(num_cons), lb=lb, ub=ub)
        end
    end
    return prob
end

"""
---
HyperelasticProblem(data::Vector{HyperelasticData}, model, u₀, ps; loss=L2DistLoss(), agg=AggMode.Mean(), cons=(x, p) -> [true], kwargs...)
---
Returns an `OptimizationProblem` for solving with GalacticOptim.jl. `data` is a vector of hyperelastic data and fits to all sets equally, `model` is the strain energy density as a function of the parameters (i.e. `f(p) = W(p)(λ⃗)`). `ps` is any hyperparameters for the model (currently not supported). `loss` defines the loss function to be used in the optimization. Currently defaults to the ``L^2``-norm between the predicted and experimental data. `agg` defines the aggregration mode of the errors, defaults to the mean of the errors. `cons` define any constrain equations involving the parameters of the model. `kwargs` are passed to `OptimizationProblem`. To set parameter bounds, use the keywords `lb` and `ub` respectively.

"""
function HyperelasticProblem(data::Vector{AbstractHyperelasticData}, model, u₀, ps; loss=L2DistLoss(), agg=AggMode.Mean(), cons=(x, p) -> [true], kwargs...)

    s = getindex.(vcat(collect.(zip(getfield.(data, :s⃗)...))...), 1) |> transpose
    λ = collect.(vcat(collect.(zip(getfield.(data, :λ⃗)...))...))

    stresses_provided = size(s, 1)

    function ŝ(p)
        s⃗ = map(x -> SecondPiolaKirchoffStressTensor(ψ, x, p), collect.(data.λ⃗))
        Δs = [s⃗[1] - s⃗[3], s⃗[2] - s⃗[3], s⃗[1] - s⃗[2]]
        return hcat(Δs...)[1:stresses_provided, :]
    end

    f(p, _) = [value(loss, s, ŝ(p), agg)]

    cons = constraints(ψ, data)
    lb, ub = parameter_bounds(ψ, data)
    func = OptimizationFunction(f, Optimization.AutoForwardDiff())
    prob = OptimizationProblem(func, u₀, ps; kwargs...)
    if !isnothing(cons)
        num_cons = length(cons(u₀, ps))
        func = remake(func, cons=cons)
        prob = remake(prob, f=func, lcons=zeros(num_cons))
    elseif !isnothing(lb) || !isnothing(ub)
        prob = remake(prob, lb=lb, ub=ub)
    end
    return prob
end
