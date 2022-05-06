using Hyperelastics
using GalacticOptim
using Plots
using Optim
using CSV
using DataFrames
using LabelledArrays
using ComponentArrays

##
pgfplotsx()

## Load Data
data_file = "example/Treloar/uniaxial.csv"
df = CSV.read(data_file, DataFrame)
data = uniaxial_data(df.stress_mpa .* 1e6, df.stretch)

## Models to fit:
models = [
    ArrudaBoyce,
    EdwardVilgis,
    Gent,
    LandisVandel,
    MooneyRivlin,
    NeoHookean,
    Ogden,
]

## Inital Guesses
u₀ = Dict(
    Gent => ComponentVector(μ=200e3, Jₘ=100.0),
    ArrudaBoyce => ComponentVector(μ=243221.2739932059, N=18.29074828993436),
    MooneyRivlin => ComponentVector(C10=200e3, C01=100e3),
    Ogden => ComponentVector(μ=[20e3, 30e3, 40e3], α=[1.0, 2.0, 3.0]),
    LandisVandel => ComponentVector(μ=247e3),
    EdwardVilgis => ComponentVector(Ns = 200e3, Nc = 200e3, α = 0.1, η = 0.2),
    NeoHookean => ComponentVector(μ = 0.246e6)
)

## Solution object
sol = Dict{Function,ComponentVector}()

## Fit the models
p = scatter(getindex.(data.λ⃗, 1), getindex.(data.s⃗, 1) ./ 1e6, label="Experimental", xlims=(1, 8), ylims = (0, 6), xlabel="Stretch (λ)", ylabel="Engineering Stress [MPa]", legend=:topleft)

for model in models
    display(model)
    if model == Gent
        Jₘ_min = sum(collect.(data.λ⃗)[end] .^ (2)) - 3
        lb = ComponentVector(μ=0.0, Jₘ=Jₘ_min)
        ub = ComponentVector(μ=Inf, Jₘ=Inf)
        HEProblem = HyperelasticProblem(data, model, u₀[model], [], lb=lb, ub=ub)
        sol[model] = solve(HEProblem, LBFGS()).u
    elseif model == MooneyRivlin
        lb = ComponentVector(C10=0.0, C01=0.0)
        ub = ComponentVector(C10=Inf, C01=Inf)
        HEProblem = HyperelasticProblem(data, model, u₀[model], [], lb=lb, ub=ub)
        sol[model] = solve(HEProblem, BFGS()).u
    elseif model == EdwardVilgis
        I₁_max = maximum(I₁.(collect(data.λ⃗)))
        α_max = sqrt(1/I₁_max)
        lb = ComponentVector(Ns=0, Nc=0, α=00.0, η=0.0)
        ub = ComponentVector(Ns=Inf, Nc=Inf, α=α_max, η=Inf)
        HEProblem = HyperelasticProblem(data, model, u₀[model], [], lb=lb, ub=ub)
        sol[model] = solve(HEProblem, LBFGS()).u
    else
        HEProblem = HyperelasticProblem(data, model, u₀[model], [])
        sol[model] = solve(HEProblem, LBFGS()).u
    end
    s⃗ = s⃗̂(model, sol[model], collect.(data.λ⃗))
    s₁ = getindex.(s⃗, 1)
    plot!(p, getindex.(data.λ⃗, 1), s₁ ./ 1e6, label=string(model))
end
display(sol)
p
savefig(p, "TreloarData.tex")