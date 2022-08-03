"""
---
HyperelasticProblem(data::HyperelasticData, model, u₀, ps; loss=L2DistLoss(), agg=AggMode.Mean(), cons=(x, p) -> [true], kwargs...)
---
Returns an `OptimizationProblem` for solving with GalacticOptim.jl. `data` is the hyperelastic experimental data, `model` is the strain energy density as a function of the parameters (i.e. `f(p) = W(p)(λ⃗)`). `ps` is any hyperparameters for the model (currently not supported). `loss` defines the loss function to be used in the optimization. Currently defaults to the ``L^2``-norm between the predicted and experimental data. `agg` defines the aggregration mode of the errors, defaults to the mean of the errors. `cons` define any constrain equations involving the parameters of the model. `kwargs` are passed to `OptimizationProblem`. To set parameter bounds, use the keywords `lb` and `ub` respectively.

"""
function HyperelasticProblem(data::AbstractHyperelasticData, ψ::AbstractHyperelasticModel, u₀, ps; loss=L2DistLoss(), agg=AggMode.Mean(), ad = Optimization.AutoForwardDiff(), kwargs...)
    s = hcat(collect.(data.s⃗)...)

    stresses_provided = size(s, 1)

    function ŝ(p)
        s⃗ = NominalStressFunction(ψ, p)
        s₁₂₃ = @. s⃗(collect(data.λ⃗))
        return hcat(s₁₂₃...)[1:stresses_provided, :]
    end

    f(p, _) = [value(loss, s, ŝ(p), agg)]

    cons = constraints(ψ, data)
    lb, ub = parameter_bounds(ψ, data)
    @show lb, ub
    func = OptimizationFunction(f, ad)
    prob = OptimizationProblem(func, u₀, ps)
    # Check for Constraints
    if !isnothing(cons)
        println("Has Constraints")
        num_cons = length(cons(u₀, ps))
        @show num_cons
        @show cons(u₀, ps)
        func = OptimizationFunction(f, ad, cons=cons)
        prob = OptimizationProblem(func, u₀, ps, lcons = zeros(num_cons))
    end
    # Check for Bounds
    if !isnothing(lb) || !isnothing(ub)
        println("Has Bounds")
        ax = Axis(Hyperelastics.parameters(ψ))
        if !isnothing(lb)&&!isnothing(ub)
            lb = ComponentVector(lb)
            ub = ComponentVector(ub)
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
        s⃗ = NominalStressFunction(ψ, p)
        s₁₂₃ = @. s⃗(collect(λ))
        return hcat(s₁₂₃...)[1:stresses_provided, :]
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
