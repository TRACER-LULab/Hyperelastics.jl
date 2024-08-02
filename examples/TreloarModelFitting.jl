using Hyperelastics
using Optimization, OptimizationOptimJL
using ComponentArrays: ComponentVector
using DifferentiationInterface
import ForwardDiff
using CairoMakie, MakiePublication
set_theme!(theme_latexfonts())
f = Figure(size=(800, 800))
ax = Makie.Axis(f[1, 1], xlabel="Stretch [-]", ylabel="Stress [kg/cm²]")
treloar_data = Treloar1944Uniaxial()
scatter!(ax,
    getindex.(treloar_data.data.λ, 1),
    getindex.(treloar_data.data.s, 1),
    label="Treloar 1944 Experimental",
    color=:black
)
axislegend(position=:lt)
save("treloar_data.png", f)
models = Dict(
    Gent => ComponentVector(μ=240e-3, J_m=80.0),
    EdwardVilgis => ComponentVector(Ns=0.10, Nc=0.20, α=0.001, η=0.001),
    NeoHookean => ComponentVector(μ=200e-3),
    Beda => ComponentVector(C1=0.1237, C2=0.0424, C3=7.84e-5, K1=0.0168, α=0.9, β=0.68, ζ=3.015)
)

sol = Dict{Any,SciMLSolution}()
for (ψ, p_0) in models
    HEProblem = HyperelasticProblem(ψ(), treloar_data, p_0, ad_type=AutoForwardDiff())
    sol[ψ] = solve(HEProblem, NelderMead())
end
return sol # hide
f = Figure(size=(400, 400))
ax = Makie.Axis(f[1, 1], xlabel="Stretch [-]", ylabel="Stress [kg/cm²]")
for (ψ, p) in sol
    pred = predict(ψ(), treloar_data, p.u, ad_type=AutoForwardDiff())
    lines!(ax, getindex.(pred.data.λ, 1), getindex.(pred.data.s, 1), label=string(ψ))
end
scatter!(ax, getindex.(treloar_data.data.λ, 1), getindex.(treloar_data.data.s, 1), label="Treloar 1944 Experimental", color=:black)
axislegend(position=:lt)
f
save("treloar_data_fits.png", f) # hide
