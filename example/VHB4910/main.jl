using Hyperelastics
using GalacticOptim
using Optim
using DataFrames
using ForwardDiff
using CSV
using LabelledArrays
using Plots

data_file = "example/VHB4910/monotonic_test_material_Heated-VHB4910_specimen_Dumbbell_H_strainrate_0.3_length_25.0_width_4.55_thickness_1.0.csv"
data_file = "example/VHB4910/UCN Longitudinal Uniaxial Test Data.csv"
data_file = "example/VHB4910/UCN Transverse Uniaxial Test Data.csv"
##
df = CSV.read(data_file, DataFrame)
data = uniaxial_data(df.Stress*1e6, df.Strain.+1)
u₀ = LVector(μ = 24e3)
HEProblem = HyperelasticProblem(data, NeoHookean, u₀, [])
sol = solve(HEProblem, NelderMead())
W = NeoHookean(sol.u)
∂s = ForwardDiff.gradient.(W, collect.(data.λ⃗))
λ₁ = getindex.(collect.(data.λ⃗), 1)
λ₃ = getindex.(collect.(data.λ⃗), 3)
s₁ = getindex.(collect.(∂s), 1)
s₃ = getindex.(collect.(∂s), 3)
σ11 = λ₁.*s₁
σ33 = λ₃.*s₃
s₁ = (σ11.-σ33)./λ₁

plot(getindex.(data.λ⃗, 1), s₁)
plot!(getindex.(data.λ⃗, 1),getindex.(data.s⃗, 1))