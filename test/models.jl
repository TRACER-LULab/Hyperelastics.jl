@testset failfast = true showtiming = true "Hyperelastics Models" begin
    # AD Backends to test.
    ADs = [AutoForwardDiff(), AutoFiniteDiff(), AutoZygote(), AutoEnzyme()]

    # Determine if the model is exported by hyperelastics.
    usemodel(model) = Base.isexported(Hyperelastics, Symbol(model))

    # collect all incompressible hyperelastic models
    incompressible_models = filter(usemodel, subtypes(Hyperelastics.AbstractIncompressibleModel))

    # Collect all compressible hyperelastics models
    compressible_models = filter(usemodel, subtypes(Hyperelastics.AbstractCompressibleModel))

    # Collect all available incompressible hyperelastic models with invariant forms
    invariant_incompressible_models = filter(Base.Fix2(applicable, InvariantForm()), incompressible_models)

    # Test the incompressible form of the model
    for model in incompressible_models
        # Instantiate model
        ψ = model()
        @show model
        @test ψ isa Hyperelastics.AbstractIncompressibleModel

        # Create an empty parameter set
        guess = Dict{Symbol,Union{Matrix{Float64},Vector{Float64},Float64}}()

        # Get the parametrs for the model and check return type
        ps = parameters(ψ)
        @test ps isa Tuple

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
                guess[k] = !isinf(v) && guess[k] < v ? (float(v) + 0.9) : (1.0)
            end
        elseif !isnothing(ub)
            for (k, v) in pairs(ub)
                guess[k] = !isinf(v) && guess[k] > v ? (float(v) - 0.9) : (1.0)
            end
        end
        guess = NamedTuple(guess)
        # Example deformation for testing
        λ⃗ = [1.1, inv(sqrt(1.1)), inv(sqrt(1.1))]
        F = diagm(λ⃗)
        λ⃗_c = λ⃗ ./ 0.99
        F_c = diagm(λ⃗_c)

        for compressible_model in compressible_models
            ψ̄ = compressible_model(ψ)
            @test ψ̄ isa Hyperelastics.AbstractCompressibleModel

            compressible_guess = (κ=1.1, ψ=guess)

            for compressible_deformation in [λ⃗_c, F_c]

                # Principal Value Form Test
                W = StrainEnergyDensity(ψ̄, compressible_deformation, compressible_guess)
                @test !isnan(W)
                @test !isinf(W)

                for AD in ADs
                    # Second Piola Kirchoff Test
                    s = SecondPiolaKirchoffStressTensor(ψ̄, compressible_deformation, compressible_guess; ad_type=AD)
                    @test sum(isnan.(s)) == 0
                    @test sum(isinf.(s)) == 0

                    # Cauchy Stress Test
                    σ = CauchyStressTensor(ψ̄, compressible_deformation, compressible_guess; ad_type=AD)
                    @test sum(isnan.(σ)) == 0
                    @test sum(isinf.(σ)) == 0

                    # Predict Test
                    @test predict(ψ̄, Treloar1944Uniaxial(), compressible_guess, ad_type=AD) isa Hyperelastics.HyperelasticUniaxialTest
                    @test predict(ψ̄, Kawabata1981(1.04), compressible_guess, ad_type=AD) isa Hyperelastics.HyperelasticBiaxialTest
                    # Strain Energy Density test for deformation gradient matrix
                    if model in invariant_incompressible_models && compressible_deformation isa Matrix
                        ψ̄_inv = compressible_model(model(InvariantForm()))

                        # Strain Energy Density Test
                        W = StrainEnergyDensity(ψ̄_inv, compressible_deformation, compressible_guess)
                        @test !isnan(W)
                        @test !isinf(W)

                        # Second Piola Kirchoff Test
                        s = SecondPiolaKirchoffStressTensor(ψ̄_inv, compressible_deformation, compressible_guess; ad_type=AD)
                        @test sum(isnan.(s)) == 0
                        @test sum(isinf.(s)) == 0

                        # Cauchy Stress Test
                        σ = CauchyStressTensor(ψ̄_inv, compressible_deformation, compressible_guess; ad_type=AD)
                        @test sum(isnan.(σ)) == 0
                        @test sum(isinf.(σ)) == 0
                    end
                end

                # Invariant Form Test
                if model in invariant_incompressible_models
                    ψ̄_inv = compressible_model(model(InvariantForm()))
                    # Strain Energy Density for vector of invariants
                    W = StrainEnergyDensity(ψ̄_inv, [I₁(compressible_deformation), I₂(compressible_deformation), I₃(compressible_deformation)], compressible_guess)
                    @test !isnan(W)
                    @test !isinf(W)
                end
            end
        end
    end

    # Test for unimplemented functions
    struct TestModel{T} <: Hyperelastics.AbstractIncompressibleModel{T}
        TestModel(::R=PrincipalValueForm()) where {R} = new{R}()
    end
    ψ = TestModel()
    ψ_inv = TestModel(InvariantForm())
    @test_throws ArgumentError parameters(ψ)
    @test_throws ArgumentError StrainEnergyDensity(ψ, ones(3), ())
    @test_throws ArgumentError StrainEnergyDensity(ψ_inv, ones(3), ())

    @test_throws ArgumentError Hyperelastics.∂ψ(ψ, ones(3), (), nothing)
end
