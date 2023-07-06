module HyperelasticsUnitfulExt

using Unitful
using Optimization
using ComponentArrays
using Hyperelastics
using LossFunctions
using Statistics
using ADTypes
using Optimization.SciMLBase

display("Due to the state of the Unitful ecosystem, the only differentiation package supported by Hyperelastics.jl with Unitful.jl is FiniteDiff.jl. Use `AutoFiniteDiff()` to use the differentiation methods from FiniteDiff.jl")

function Hyperelastics.HyperelasticProblem(
    ψ::Hyperelastics.AbstractHyperelasticModel,
    test::Hyperelastics.AbstractHyperelasticTest{T,S},
    u0::U,
    adtype::AutoFiniteDiff,
    loss=L2DistLoss();
    lb=parameter_bounds(ψ, test).lb,
    ub=parameter_bounds(ψ, test).ub,
    int=nothing,
    lcons=nothing,
    ucons=nothing,
    sense=nothing,
    kwargs...) where {T,S<:Quantity,U}

    display("UNIT VERSION")
    a = typeof(adtype)
    b = Symbol(split(split(string(a), "{")[1], ".")[end])

    optimization_AD_type = getfield(Optimization, b)()

    function f(ps, (; ψ, test, kwargs, loss, units))
        pred = predict(ψ, test, ps.*units, adtype, kwargs...)
        # res = sum(abs, value(loss, test.data.s, pred.data.s, AggMode.Mean()))
        res = map(i -> L2DistLoss().(i[1], i[2]), zip(pred.data.s, test.data.s)) |> mean

        res = ustrip.(res)
        return res
    end

    lb, ub = parameter_bounds(ψ, test)

    if !isnothing(lb) && !isnothing(ub)
        lb = ComponentVector(lb)
        ub = ComponentVector(ub)
    elseif !isnothing(lb)
        lb = ComponentVector(lb)
        ub = u0 .* Inf
    elseif !isnothing(ub)
        ub = ComponentVector(ub)
        lb = u0 .* -Inf
    else
        ub = u0 .* Inf
        lb = u0 .* -Inf
    end

    model_ps = parameters(ψ)

    for p in model_ps
        if !isnothing(lb)
            if (u0[p] < lb[p])
                @error "Parameter $p = $(u₀[p]) is less than lower bound of $(lb[p])"
                return nothing
            end
        end
        if !isnothing(ub)
            if (u0[p] > ub[p])
                @error "Parameter $p = $(u₀[p]) is greater than upper bound of $(ub[p])"
                return nothing
            end
        end
    end

    func = OptimizationFunction(f, optimization_AD_type)
     # Check for Bounds
    p = (ψ=ψ, test=test, loss=loss, kwargs=kwargs, units = unit.(u0))
    u0 = ustrip.(u0)
    lb = (isnothing(lb)) ? (lb) : (ustrip.(lb))
    ub = (isnothing(ub)) ? (ub) : (ustrip.(ub))
    HyperelasticProblem{isinplace(func),typeof(func),typeof(u0),typeof(p),
        typeof(lb),typeof(ub),typeof(int),typeof(lcons),typeof(ucons),
        typeof(sense),typeof(kwargs)}(func, u0, p, lb, ub, int, lcons, ucons, sense,
        kwargs)
end

function Optimization.SciMLBase.solve(
    prob::HyperelasticProblem,
    alg,
    args...;
    kwargs...)::Optimization.SciMLBase.AbstractOptimizationSolution

    conv_prob = OptimizationProblem(prob.f, prob.u0, prob.p, lb=prob.lb, ub=prob.ub, int=prob.int, lcons=prob.lcons, ucons=prob.ucons, sense=prob.sense, prob.kwargs...)

    sol = Optimization.SciMLBase.solve(conv_prob, alg, args...; kwargs...)
    Optimization.SciMLBase.OptimizationSolution{
        typeof(sol.retcode),
        ndims(sol.u),
        typeof(sol.u),
        typeof(sol.cache),
        typeof(sol.alg),
        typeof(sol.objective),
        typeof(sol.original),
        typeof(sol.solve_time),
        typeof(sol.stats)}(
        sol.u,
        sol.cache,
        sol.alg,
        sol.objective,
        sol.retcode,
        sol.original,
        sol.solve_time,
        sol.stats)
end

end
