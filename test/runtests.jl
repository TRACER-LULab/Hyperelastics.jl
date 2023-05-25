using LinearAlgebra
using Hyperelastics
using Test
using ForwardDiff
using FiniteDiff
using InteractiveUtils

@testset "Hyperelastics.jl" begin
    # Test loading the prommvided Datasets
    treloar_data = Treloar1944Uniaxial()
    # Skip Data-Driven and Compressible Models
    compressible_models = [
        GeneralCompressible,
        LogarithmicCompressible,
    ]
    # Check to see if it is a publicly exported model and that it is not one of the compresible model wrappers.
    usemodel(model) = !Base.Fix2(in, compressible_models)(model) && Base.isexported(Hyperelastics, Symbol(model))
    models = filter(usemodel, subtypes(Hyperelastics.AbstractHyperelasticModel))
    # Get all available analytical hyperelastic models
    for model in models
        # Check Initialization of model for default settings
        ψ = model()
        # Check that parameters are defined for model
        ps = parameters(ψ)
        # Check for bounds on parameters
        lb, ub = parameter_bounds(ψ, treloar_data)
        # Setup a set of parameters to check with
        guess = Dict{Symbol, Union{Vector{Float64},Float64}}()
        for p in ps
            if contains(string(p), '⃗')
                guess[p] = ones(10)
            else
                guess[p] = 1.0
            end
        end
        if !isnothing(lb) && !isnothing(ub)
            for (k,v) in pairs(lb)
                lb_val = !isinf(getfield(lb, k)) ? (float(getfield(lb, k))) : (1.0)
                ub_val = !isinf(getfield(ub, k)) ? (float(getfield(ub, k))) : (1.0)
                guess[k] = (lb_val + ub_val) / 2.0
            end
        elseif !isnothing(lb)
            for (k,v) in pairs(lb)
                guess[k] = !isinf(getfield(lb, k)) ? (getfield(lb, k))+1.0 : (1.0)
            end
        elseif !isnothing(ub)
            for (k, v) in pairs(ub)
                guess[k] = !isinf(getfield(ub, k)) ? (getfield(ub, k))-1.0 : (1.0)
            end
        end

        # Create the parameter set to test against
        guess = NamedTuple(guess)

        # Test deformations.
        λ⃗ = [1.1, inv(sqrt(1.1)), inv(sqrt(1.1))]
        F = diagm(λ⃗)
        λ⃗_c = copy(λ⃗)
        λ⃗_c[2] = λ⃗_c[3] = λ⃗_c[2]*sqrt(0.97)
        F_c = diagm(λ⃗_c)
        # Test the incompressible form of the model
        for deformation in [λ⃗, F]
            @test !isnan(StrainEnergyDensity(ψ, deformation, guess))
            @test !isinf(StrainEnergyDensity(ψ, deformation, guess))
        end
        # Test against all compressible forms of the model
        for compressible_model in compressible_models
            for compressible_deformation in [λ⃗_c, F_c]
                compressible_guess = (κ = 1.1, ψ = guess)
                # Create the compressible model
                ψ̄ = compressible_model(ψ)
                # Strain Energy Density Test
                @test !isnan(StrainEnergyDensity(ψ̄, compressible_deformation, compressible_guess))
                @test !isinf(StrainEnergyDensity(ψ̄, compressible_deformation, compressible_guess))
                # Loop over AD Backends
                # for AD in [AutoForwardDiff()]
                #     # Second Piola Kirchoff Stress Test
                #     SecondPiolaKirchoffStressTensor(ψ̄, deformation, compressible_guess, AD)
                #     # Cauchy Stress Tensor Test
                # end
            end

        end
    end
end
