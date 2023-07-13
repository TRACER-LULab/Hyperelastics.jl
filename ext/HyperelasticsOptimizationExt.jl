module HyperelasticsOptimizationExt

using Optimization
using ComponentArrays
using Hyperelastics
using LossFunctions
using Statistics
using ADTypes
# using Optimization.SciMLBase

export HyperelasticProblem

# include(joinpath(@__DIR__, "..", "src", "macro_micro_macro_model.jl"))

# Optimization.SciMLBase.@add_kwonly function HyperelasticProblem{iip}(
#     f::OptimizationFunction{iip},
#     u0,
#     p=NullParameters();
#     lb=nothing, ub=nothing, int=nothing,
#     lcons=nothing, ucons=nothing,
#     sense=nothing, kwargs...) where {iip}

#     if xor(lb === nothing, ub === nothing)
#         error("If any of `lb` or `ub` is provided, both must be provided.")
#     end

#     OptimizationProblem(f, u0, p, lb, ub, int, lcons, ucons, sense,
#             kwargs)
# end

function Hyperelastics.HyperelasticProblem(
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

    function f(ps, p)
        ψ, test, loss, ad_type, kwargs = p
        pred = predict(ψ, test, ps; ad_type, kwargs...)
        res = map(i -> loss.(i[1], i[2]), zip(pred.data.s, test.data.s)) |> mean
        return res
    end

    u0 = ComponentVector(u0)
    if !isnothing(lb) && !isnothing(ub)
        lb = ComponentVector(lb,)
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
    p = (ψ, test, loss, ad_type, kwargs)
    OptimizationProblem(
        func,
        u0,
        p;
        lb,
        ub,
        int,
        lcons,
        ucons,
        sense,
    )
end

function Hyperelastics.HyperelasticProblem(
    ψ::Hyperelastics.AbstractHyperelasticModel,
    tests::Vector{R},
    u0;
    ad_type::ADTypes.AbstractADType,
    loss=L2DistLoss(),
    lb=parameter_bounds(ψ, tests).lb,
    ub=parameter_bounds(ψ, tests).ub,
    int=nothing,
    lcons=nothing,
    ucons=nothing,
    sense=nothing,
    kwargs...) where {R<:Hyperelastics.AbstractHyperelasticTest}

    get_s(test) = hcat(test.data.s...)

    function f(ps, p)
        ψ, tests, loss, kwargs, ad_type = p
        preds = predict(ψ, tests, ps; ad_type, kwargs...)

        s = get_s.(tests)
        ŝ = get_s.(preds)
        res = map(idx -> mean(map(i-> loss.(i[1], i[2]), zip(ŝ[idx], s[idx]))), eachindex(s))|>mean
        return res
    end

    u0 = ComponentVector(u0)
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
    p = (ψ, tests, loss, kwargs, ad_type)

    OptimizationProblem(
        func,
        u0,
        p;
        lb,
        ub,
        int,
        lcons,
        ucons,
        sense,
    )
end

# function Optimization.SciMLBase.solve(
#     prob::HyperelasticProblem,
#     alg,
#     args...;
#     kwargs...)::Optimization.SciMLBase.AbstractOptimizationSolution

#     conv_prob = OptimizationProblem(prob.f, prob.u0, prob.p, lb=prob.lb, ub=prob.ub, int=prob.int, lcons=prob.lcons, ucons=prob.ucons, sense=prob.sense, prob.kwargs...)

#     Optimization.SciMLBase.solve(conv_prob <: AbstractOptimizationProblem, alg, args...; kwargs...)::Optimization.SciMLBase.AbstractOptimizationSolution
# end

end
