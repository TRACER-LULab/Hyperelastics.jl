@testset "Model Fitting" begin
    using Optimization, OptimizationOptimJL, ComponentArrays, ForwardDiff, NaNMath
    # Determine if the model is exported by hyperelastics.
    usemodel(model) = Base.isexported(Hyperelastics, Symbol(model))

    # collect all incompressible hyperelastic models
    incompressible_models =
        filter(usemodel, subtypes(Hyperelastics.AbstractIncompressibleModel))

    # Collect all compressible hyperelastics models
    compressible_models =
        filter(usemodel, subtypes(Hyperelastics.AbstractCompressibleModel))

    # Collect all available incompressible hyperelastic models with invariant forms
    invariant_incompressible_models =
        filter(Base.Fix2(applicable, InvariantForm()), incompressible_models)

    # Create initial


    # Test the incompressible form of the model
    for model in [NeoHookean]
        # Test the incompressible form of the model
        ψ = model()
        @show model
        # Get the parameters for the model and check return type
        ps = parameters(ψ)
        # Create an empty parameter set
        guess = Dict{Symbol,Union{Matrix{Float64},Vector{Float64},Float64}}()
        # Use scalars, vectors, or matrices based on the symbol provided.
        for p in ps
            if ψ isa GeneralMooneyRivlin
                guess[p] = [1.0 1.0; 1.0 1.0]
            elseif ψ isa HorganMurphy
                guess[p] = if p == :μ
                    1.0
                elseif p == :c
                    2.0
                elseif p == :Jₘ
                    200.0
                end
            elseif ψ isa VanDerWaals
                guess[p] = if p == :μ
                    0.355
                elseif p == :λm
                    10.0
                elseif p == :β
                    1.0
                elseif p == :α
                    0.25
                end
            elseif contains(string(p), '⃗')
                guess[p] = ones(10)
            else
                guess[p] = 1.0
            end
        end

        # Check for bounds on the model
        lb, ub = parameter_bounds(ψ, Treloar1944Uniaxial())
        @test lb isa NamedTuple || isnothing(lb)
        @test ub isa NamedTuple || isnothing(ub)

        lb, ub = parameter_bounds(ψ, [Treloar1944Uniaxial(), Kawabata1981(1.04)])
        @test lb isa NamedTuple || isnothing(lb)
        @test ub isa NamedTuple || isnothing(ub)

        # Move the guess to within the parameter bounds
        if ψ isa VanDerWaals || ψ isa HorganMurphy
            nothing
        elseif !isnothing(lb) && !isnothing(ub)
            for (k, v) in pairs(lb)
                lb_val = !isinf(v) && guess[k] < v ? (float(v) + 0.9) : (1.0)
                ub_val = !isinf(v) && guess[k] > v ? (float(v) - 0.9) : (1.0)
                guess[k] = (lb_val + ub_val) / 2.0
            end
        elseif !isnothing(lb)
            for (k, v) in pairs(lb)
                @show v
                guess[k] = !isinf(v) && guess[k] < v ? (float(v) + 0.9) : (1.0)
            end
        elseif !isnothing(ub)
            for (k, v) in pairs(ub)
                guess[k] = !isinf(v) && guess[k] > v ? (float(v) - 0.9) : (1.0)
            end
        end
        guess = ComponentVector(guess)
        display(guess)
        prob =
            HyperelasticProblem(ψ, Treloar1944Uniaxial(), guess; ad_type = AutoFiniteDiff())
        display(prob.lb)
        sol = solve(prob, LBFGS())
        @test sol.retcode == ReturnCode.Success

        # sol = solve(prob, LBFGS())
        # @test sol.retcode == ReturnCode.Success
    end
end
