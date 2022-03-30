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
##
df = CSV.read(data_file, DataFrame)
data = uniaxial_data(df.Stress*1e6, df.Strain.+1)
u₀ = LVector(μ = 45e6, Jₘ = 0.25)
HEProblem = HyperelasticProblem(data, Gent, u₀, [])
sol = solve(HEProblem, NelderMead())
W = Gent(sol.u)
s⃗ = ForwardDiff.gradient.(W, collect.(data.λ⃗))
s₁ = getindex.(s⃗, 1).-getindex.(s⃗, 3)
plot(getindex.(data.λ⃗, 1), s₁)
plot!(getindex.(data.λ⃗, 1),getindex.(data.s⃗, 1))