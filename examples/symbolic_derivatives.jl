using Hyperelastics
using Symbolics
using Base.Threads
using ProgressMeter
##
# Get Models
required = [
    :NominalStressFunction,
    :TrueStressFunction,
    :parameters,
    :StrainEnergyDensityFunction,
    :AbstractHyperelasticData,
    :BiaxialHyperelasticData,
    :UniaxialHyperelasticData,
    :HyperelasticProblem,
    :Hyperelastics,
    :I₁,
    :I₂,
    :I₃,
    :J,
    :Kawabata1981
]

data_driven = [
    :MacroMicroMacro,
    :SussmanBathe,
    :DataDrivenAverageChainBehavior,
    :SussmanBathe,
    :StabilitySmoothedSussmanBathe,
    :LatoreeMontans
]

hyperelastic_models = filter(x -> !(x in required) && !(x in data_driven), names(Hyperelastics))
hyperelastic_models = map(x -> getfield(Hyperelastics, Symbol(x)), hyperelastic_models)
##
@variables λ₁, λ₂, λ₃
λ⃗ = [λ₁, λ₂, λ₃]
# @variables λ⃗[1:3]
D = Differential.(λ⃗)
# D = Differential.(Symbolics.scalarize(λ⃗))
# open("generated_symbolic_derivatives.jl", "w") do f
# end
# open("generated_symbolic_derivatives.jl", "a") do f
@threads for model in hyperelastic_models
    @show model
    if isfile(joinpath(@__DIR__, "..", "src", "GeneratedNominalStressFunctions", string(typeof(model())) * ".jl"))
        display("skipping $(model())")
        continue
    end
    if hasfield(model, :ℒinv)
        display("skipping $(model())")
        continue
    end
    if model() == GeneralMooneyRivlin()
        display("skipping $(model())")
        continue
    end
    symbol_ps = Hyperelastics.parameters(model())
    if sum(contains.(string.(symbol_ps), "⃗")) > 0
        display("skipping $(model())")
        continue
    end
    ps = symbol_ps .|> Symbolics.Sym{Number}
    p = Dict(symbol_ps .=> ps) |> NamedTuple
    W = StrainEnergyDensityFunction(model(), p)
    ∂W∂λ = [D[1](W(Symbolics.scalarize(λ⃗))), D[2](W(Symbolics.scalarize(λ⃗))), D[3](W(Symbolics.scalarize(λ⃗)))] .|> expand_derivatives
    s⃗ = ∂W∂λ .- ∂W∂λ[3] .* λ⃗[3] ./ λ⃗#  .|> simplify
    σ⃗ = s⃗ .* λ⃗ #.|> simplify
    s⃗_expr = build_function(
        s⃗,
        [λ⃗[1], λ⃗[2], λ⃗[3], ps...],
        expression=Val{true}
    )
    σ⃗_expr = build_function(
        σ⃗,
        [λ⃗[1], λ⃗[2], λ⃗[3], ps...],
        linenumbers=false,
        expression=Val{true}
    )
    parameter_string1 = string(ps)
    parameter_string1 = parameter_string1[1] * ";" * parameter_string1[2:end]
    parameter_string2 = string(ps)
    parameter_string2 = replace(parameter_string2, "(" => "", ")" => "")
    open(joinpath(@__DIR__, "..", "src", "GeneratedNominalStressFunctions", string(typeof(model())) * ".jl"), "w") do f
        write(f, "function NominalStressFunction(ψ::$(typeof(model())), $(parameter_string1))\n")
        write(f, "func = ")
        write(f, string(s⃗_expr[1]))
        write(f, "\ns⃗(λ⃗) = func([λ⃗[1], λ⃗[2], λ⃗[3], $(parameter_string2)])")
        write(f, "\nend\n")
    end
    println("---\tWriting σ function\t ---")
    open(joinpath(@__DIR__, "..", "src", "GeneratedTrueStressFunctions", string(typeof(model())) * ".jl"), "w") do f
        write(f, "function TrueStressFunction(ψ::$(typeof(model())), $(parameter_string1))\n")
        write(f, "func = ")
        write(f, string(σ⃗_expr[1]))
        write(f, "\nσ⃗(λ⃗) = func([λ⃗[1], λ⃗[2], λ⃗[3], $(parameter_string2)])")
        write(f, "\nend\n")
    end
end
