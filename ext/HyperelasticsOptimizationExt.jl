module HyperelasticsOptimizationExt

using Optimization
using ComponentArrays
using Hyperelastics
using LossFunctions
using Statistics
using ADTypes
using Optimization.SciMLBase

export HyperelasticProblem

Optimization.SciMLBase.@add_kwonly function HyperelasticProblem{iip}(
    f::OptimizationFunction{iip},
    u0,
    p=NullParameters();
    lb=nothing, ub=nothing, int=nothing,
    lcons=nothing, ucons=nothing,
    sense=nothing, kwargs...) where {iip}

    if xor(lb === nothing, ub === nothing)
        error("If any of `lb` or `ub` is provided, both must be provided.")
    end

    HyperelasticProblem{iip,typeof(f),typeof(u0),typeof(p),
            typeof(lb),typeof(ub),typeof(int),typeof(lcons),typeof(ucons),
            typeof(sense),typeof(kwargs)}(f, u0, p, lb, ub, int, lcons, ucons, sense,
            kwargs)
end

function HyperelasticProblem(
    ψ::Hyperelastics.AbstractHyperelasticModel,
    test::Hyperelastics.AbstractHyperelasticTest{T,S},
    u0;
    ad_type::ADTypes.AbstractADType,
    loss=L2DistLoss(),
    lb=parameter_bounds(ψ, test).lb,
    ub=parameter_bounds(ψ, test).ub,
    int=nothing,
    lcons=nothing,
    ucons=nothing,
    sense=nothing,
    kwargs...) where {T,S}

    # optimization_AD_type = ad_type

    function f(ps, (; ψ, test, kwargs, loss))
        pred = predict(ψ, test, ps; ad_type, kwargs...)
        # @show test.data.s[1]

        # res = sum(abs, mean(loss.(pred.data.s, test.data.s)))
        res = map(i -> L2DistLoss().(i[1], i[2]), zip(pred.data.s, test.data.s)) |> mean
        return res
    end

    # lb, ub = parameter_bounds(ψ, test)

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
                @error "Parameter $p = $(u0[p]) is less than lower bound of $(lb[p])"
                return nothing
            end
        end
        if !isnothing(ub)
            if (u0[p] > ub[p])
                @error "Parameter $p = $(u0[p]) is greater than upper bound of $(ub[p])"
                return nothing
            end
        end
    end

    func = OptimizationFunction(f, ad_type)
         # Check for Bounds
    p = (ψ=ψ, test=test, loss=loss, kwargs=kwargs)

    HyperelasticProblem{isinplace(func),typeof(func),typeof(u0),typeof(p),
            typeof(lb),typeof(ub),typeof(int),typeof(lcons),typeof(ucons),
            typeof(sense),typeof(kwargs)}(func, u0, p, lb, ub, int, lcons, ucons, sense,
            kwargs)
end

    # function HyperelasticProblem(ψ::Hyperelastics.AbstractHyperelasticModel, tests::Vector{Hyperelastics.AbstractHyperelasticTest{T,S}}, u₀, adtype::ADTypes.AbstractADType, p=NullParameters(); loss=L2DistLoss(), lb=nothing, ub=nothing, int=nothing,
    #     lcons=nothing, ucons=nothing,
    #     sense=nothing, kwargs...) where {T,S}

    #     get_s(test) = hcat(test.data.s...)
    #     optimization_AD_type = getfield(Optimization, Symbol(split(string(typeof(adtype)), "{")[1]))()
    #     ps = (ψ=ψ, tests=tests, loss=loss, kwargs=kwargs)

    #     function f(ps, p)
    #         (; ψ, tests, kwargs, loss) = p
    #         preds = map(test -> predict(ψ, test, ps, adtype, kwargs...), tests)
    #         s = get_s.(tests)
    #         ŝ = get_s.(preds)
    #         res = sum(abs, map(i -> mean(value(loss, s[i], ŝ[i], AggMode.Mean())), eachindex(s)))
    #         return res
    #     end

    #     func = OptimizationFunction(f, optimization_AD_type)
    #     prob = HyperelasticProblem(func, u₀, ps)

    #     # # Check for Bounds
    #     lb, ub = parameter_bounds(ψ, tests)
    #     if !isnothing(lb) || !isnothing(ub)
    #         ax = Axis(Hyperelastics.parameters(ψ))
    #         if !isnothing(lb) && !isnothing(ub)
    #             lb = ComponentVector(lb)
    #             ub = ComponentVector(ub)
    #         elseif !isnothing(lb)
    #             lb = ComponentVector(lb)
    #             ub = ComponentVector(ones(length(lb)) * Inf, ax)
    #         elseif !isnothing(ub)
    #             ub = ComponentVector(ub)
    #             lb = ComponentVector(ones(length(ub)) * -Inf, ax)
    #         end
    #         if isnothing(prob.lcons)
    #             prob = HyperelasticProblem(func, u₀, ps, lb=lb, ub=ub)
    #         else
    #             num_cons = length(cons(u₀, ps))
    #             prob = HyperelasticProblem(func, u₀, ps, lcons=cons.lcons, ucons=cons.ucons, lb=lb, ub=ub)
    #         end
    #     end
    #     return prob
    # end

function Optimization.SciMLBase.solve(
    prob::HyperelasticProblem,
    alg,
    args...;
    kwargs...)::Optimization.SciMLBase.AbstractOptimizationSolution

    conv_prob = OptimizationProblem(prob.f, prob.u0, prob.p, lb=prob.lb, ub=prob.ub, int=prob.int, lcons=prob.lcons, ucons=prob.ucons, sense=prob.sense, prob.kwargs...)

    Optimization.SciMLBase.solve(conv_prob, alg, args...; kwargs...)::Optimization.SciMLBase.AbstractOptimizationSolution
end

end
