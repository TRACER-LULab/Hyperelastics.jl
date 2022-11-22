struct HyperelasticDataEntry{T,S}
    λ::Vector{T}
    s::Vector{S}
end

struct HyperelasticUniaxialTest <: AbstractHyperelasticTest
    data::StructVector
    name::String
    function HyperelasticUniaxialTest(λ₁, s₁; name, incompressible=true)
        @assert length(λ₁) == length(s₁) "Inputs must be the same length"
        if incompressible
            λ₂ = λ₃ = @. sqrt(1 / λ₁)
        else
            λ₂ = λ₃ = Vector{eltype(λ₁)}(undef, length(λ₁))
        end
        λ = collect.(zip(λ₁, λ₂, λ₃))
        s = collect.(zip(s₁))
        data = StructArray{HyperelasticDataEntry}((λ, s))
        new(data, name)
    end
    function HyperelasticUniaxialTest(λ₁; name, incompressible=true)
        if incompressible
            λ₂ = λ₃ = @. sqrt(1 / λ₁)
        else
            λ₂ = λ₃ = Vector{eltype(λ₁)}(undef, length(λ₁))
        end
        λ = collect.(zip(λ₁, λ₂, λ₃))
        s = collect.(zip(Vector{eltype(λ₁)}(undef, length(λ₁))))
        data = StructArray{HyperelasticDataEntry}((λ, s))
        new(data, name)
    end
end

function Base.show(io::IO, test::HyperelasticUniaxialTest)
    print(io, Term.RenderableText("Uniaxial Test: {bold} $(test.name)"))
    print(io,
        Term.Table(
            hcat(getindex.(test.data.λ, 1), getindex.(test.data.s, 1)),
            header=["λ₁", "s₁"],
            box=:ROUNDED,
        )
    )
end

struct HyperelasticBiaxialTest <: AbstractHyperelasticTest
    data::StructVector
    name::String
    function HyperelasticBiaxialTest(λ₁, λ₂, s₁, s₂; name, incompressible=true)
        @assert length(λ₁) == length(λ₂) == length(s₁) == length(s₂) "Inputs must be the same length"
        if incompressible
            λ₃ = @. 1 / λ₁ / λ₂
        else
            λ₃ = Vector{eltype(λ₁)}(undef, length(λ₁))
        end
        λ = collect.(zip(λ₁, λ₂, λ₃))
        s = collect.(zip(s₁, s₂))
        data = StructArray{HyperelasticDataEntry}((λ, s))
        new(data, name)
    end
    function HyperelasticBiaxialTest(λ₁, λ₂; name, incompressible=true)
        @assert length(λ₁) == length(λ₂) "Inputs must be the same length"
        if incompressible
            λ₃ = @. 1 / λ₁ / λ₂
        else
            λ₃ = Vector{eltype(λ₁)}(undef, length(λ₁))
        end
        s₁ = s₂ = Vector{eltype(λ₁)}(undef, length(λ₁))
        λ = collect.(zip(λ₁, λ₂, λ₃))
        s = collect.(zip(s₁, s₂))
        data = StructArray{HyperelasticDataEntry}((λ, s))
        new(data, name)
    end
end

function Base.show(io::IO, test::HyperelasticBiaxialTest)
    print(io, Term.RenderableText("Biaxial Test: {bold} $(test.name)"))
    print(io,
        Term.Table(
            hcat(getindex.(test.data.λ, 1),getindex.(test.data.λ, 2), getindex.(test.data.s, 1), getindex.(test.data.s, 2)),
            header=["λ₁", "λ₂", "s₁", "s₂"],
            box=:ROUNDED,
        )
    )
end

function NonlinearContinua.predict(ψ::AbstractHyperelasticModel, test::HyperelasticUniaxialTest, p)
    f(λ) = SecondPiolaKirchoffStressTensor(ψ, λ, p)
    λ = test.data.λ
    s = map(f, λ)
    s₁ = getindex.(s, 1)
    s₃ = getindex.(s, 3)
    λ₁ = getindex.(λ, 1)
    λ₃ = getindex.(λ, 3)
    Δs₁₃ = @. s₁ - s₃ * λ₃ / λ₁
    HyperelasticUniaxialTest(λ₁, Δs₁₃, name=test.name)
end

function NonlinearContinua.predict(ψ::AbstractHyperelasticModel, test::HyperelasticBiaxialTest, p)
    f(λ) = SecondPiolaKirchoffStressTensor(ψ, λ, p)
    λ = test.data.λ
    s = map(f, λ)
    s₁ = getindex.(s, 1)
    s₂ = getindex.(s, 2)
    s₃ = getindex.(s, 3)
    λ₁ = getindex.(λ, 1)
    λ₂ = getindex.(λ, 2)
    λ₃ = getindex.(λ, 3)
    Δs₂₃ = @. s₂ - s₃ * λ₃ / λ₁
    Δs₁₃ = @. s₁ - s₃ * λ₃ / λ₁
    HyperelasticBiaxialTest(λ₁, λ₂, Δs₁₃, Δs₂₃, name=test.name)
end

function HyperelasticProblem(ψ::AbstractHyperelasticModel, test::AbstractHyperelasticTest, u₀, ps=Nothing;
    loss=L2DistLoss(), agg=AggMode.Mean(), adtype=Optimization.AutoForwardDiff())

    function f(ps, (; ψ, test))
        ŷ = predict(ψ, test, ps)
        ŝ = getindex.(ŷ.data.s, 1)
        s = getindex.(test.data.s, 1)
        nmae = value(loss, s, ŝ, AggMode.Mean())
        return nmae
    end

    cons = constraints(ψ, test)
    lb, ub = parameter_bounds(ψ, test)
    model_ps = parameters(ψ)
    ps = (ψ=ψ, test=test)

    for p in model_ps
        if !isnothing(lb)
            (u₀[p] < lb[p]) ? (error("Parameter $p is less than lower bound of $(lb[p])")) : ()
        elseif !isnothing(ub)
            (u₀[p] > ub[p]) ? (error("Parameter $p is greater than upper bound of $(ub[p])")) : ()
        end
    end
    func = OptimizationFunction(f, adtype)
    prob = OptimizationProblem(func, u₀, ps)
    # Check for Constraints
    if !isnothing(cons)
        num_cons = length(cons.ucons)
        func = OptimizationFunction(f, adtype, cons=cons.cons)
        prob = OptimizationProblem(func, u₀, ps, lcons=cons.lcons,ucons = cons.ucons)
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
            prob = OptimizationProblem(func, u₀, ps, lcons=cons.lcons,ucons = cons.ucons, lb=lb, ub=ub)
        end
    end
    return prob
end
