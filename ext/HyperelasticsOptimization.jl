module HyperelasticsOptimization

using Optimization
using ComponentArrays
using Hyperelastics
using LossFunctions
using Statistics
using ADTypes
# export HyperelasticProblem
"""
`HyperelasticProblem(ψ::AbstractHyperelasticModel, test::AbstractHyperelasticTest, u₀, ps=Nothing;
    adb=AD.ForwardDiffBackend(), loss=L2DistLoss(), adtype=Optimization.AutoForwardDiff())`

`HyperelasticProblem(ψ::AbstractHyperelasticModel, tests::Vector{<:AbstractHyperelasticTest}, u₀, ps=Nothing;
    adb=AD.ForwardDiffBackend(), loss=L2DistLoss(), adtype=Optimization.AutoForwardDiff())`

Creates an `OptimizationProblem` for use in [`Optimization.jl`](https://docs.sciml.ai/Optimization/stable/) to find the optimal parameters.

Fields:
- `ψ`: material model to use
- `test` or `tests`: A single or vector of hyperelastics tests to use when fitting the parameters
- `u₀`: Initial guess for parameters
- `ps`: Any additional parameters for calling predict
- `adb`: Select differentiation type from [`ADTypes.jl`](https://github.com/SciML/ADTypes.jl). The type is automatically applied to the type of AD applied to the Optimization Problem also.
- `loss`: Loss function from [`LossFunctions.jl`](https://github.com/JuliaML/LossFunctions.jl)
"""
function Hyperelastics.HyperelasticProblem(ψ::Hyperelastics.AbstractHyperelasticModel, test::Hyperelastics.AbstractHyperelasticTest, u₀, adtype::ADTypes.AbstractADType; loss=L2DistLoss(), kwargs...)
    a = typeof(adtype)
    b = Symbol(split(split(string(a), "{")[1], ".")[end])
    optimization_AD_type = getfield(Optimization, b)()
    function f(ps, (; ψ, test, kwargs, loss))
        pred = predict(ψ, test, ps, adtype, kwargs...)
        s = hcat(test.data.s...)
        ŝ = hcat(pred.data.s...)
        res = sum(abs, value(loss, s, ŝ, AggMode.Mean(), ObsDim.First()))
        return res
    end

    cons = Hyperelastics.constraints(ψ, test)
    lb, ub = parameter_bounds(ψ, test)
    model_ps = parameters(ψ)
    ps = (ψ=ψ, test=test, loss=loss, kwargs=kwargs)

    for p in model_ps
        if !isnothing(lb)
            if (u₀[p] < lb[p])
                @error "Parameter $p = $(u₀[p]) is less than lower bound of $(lb[p])"
                return nothing
            end
        end
        if !isnothing(ub)
            if (u₀[p] > ub[p])
                @error "Parameter $p = $(u₀[p]) is greater than upper bound of $(ub[p])"
                return nothing
            end
        end
    end
    func = OptimizationFunction(f, optimization_AD_type)
    prob = OptimizationProblem(func, u₀, ps)
    # Check for Constraints
    if !isnothing(cons)
        num_cons = length(cons.ucons)
        func = OptimizationFunction(f, optimization_AD_type, cons=cons.cons)
        prob = OptimizationProblem(func, u₀, ps, lcons=cons.lcons, ucons=cons.ucons)
    end
    # Check for Bounds
    if !isnothing(lb) || !isnothing(ub)
        ax = Axis(Hyperelastics.parameters(ψ))
        if !isnothing(lb) && !isnothing(ub)
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
            prob = OptimizationProblem(func, u₀, ps, lcons=cons.lcons, ucons=cons.ucons, lb=lb, ub=ub)
        end
    end
    return prob
end

function Hyperelastics.HyperelasticProblem(ψ::Hyperelastics.AbstractHyperelasticModel, tests::Vector{<:Union{Hyperelastics.HyperelasticUniaxialTest,Hyperelastics.HyperelasticBiaxialTest}}, u₀, adtype::ADTypes.AbstractADType; loss=L2DistLoss(), kwargs...)

    get_s(test) = hcat(test.data.s...)
    optimization_AD_type = getfield(Optimization, Symbol(split(string(typeof(adtype)), "{")[1]))()
    ps = (ψ=ψ, tests=tests, loss=loss, kwargs=kwargs)

    function f(ps, p)
        (; ψ, tests, kwargs, loss) = p
        preds = map(test->predict(ψ, test, ps, adtype, kwargs...), tests)
        s = get_s.(tests)
        ŝ = get_s.(preds)
        res = sum(abs, map(i -> mean(value(loss, s[i], ŝ[i], AggMode.Mean(), ObsDim.First())), eachindex(s)))
        return res
    end

    func = OptimizationFunction(f, optimization_AD_type)
    prob = OptimizationProblem(func, u₀, ps)

    # # Check for Bounds
    lb, ub = parameter_bounds(ψ, tests)
    if !isnothing(lb) || !isnothing(ub)
        ax = Axis(Hyperelastics.parameters(ψ))
        if !isnothing(lb) && !isnothing(ub)
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
            prob = OptimizationProblem(func, u₀, ps, lcons=cons.lcons, ucons=cons.ucons, lb=lb, ub=ub)
        end
    end
    return prob
end
end
