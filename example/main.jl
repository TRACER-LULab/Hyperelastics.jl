using Hyperelastics
using GalacticOptim
using AbstractDifferentiation  
using ForwardDiff
using Optim
using DataFrames
using CSV
using LabelledArrays
using ComponentArrays
using PlotlyLight
using Plots

## Set the data File to read
# data_file = "example/VHB4910/monotonic_test_material_Heated-VHB4910_specimen_Dumbbell_H_strainrate_0.3_length_25.0_width_4.55_thickness_1.0.csv"
# data_file = "example/PSHBladder/UCN Longitudinal Uniaxial Test Data.csv"
# data_file = "example/PSHBladder/UCN Transverse Uniaxial Test Data.csv"
data_file = "example/Treloar/uniaxial.csv"

## Read the file and create the HyperelasticData object
df = CSV.read(data_file, DataFrame)

## Create the uniaxial test data object
# data = uniaxial_data(df.Stress*1e6, df.Strain.+1)
# data = uniaxial_data(df.stress[1:25:2550] .* 1e6, df.strain[1:25:2550] .+ 1)
data = uniaxial_data(df.stress_mpa.*1e6, df.stretch)

## GeneralizedMooneyRivlin
model = GeneralizedMooneyRivlin
u₀ = ComponentVector(C=[0.0 10e2 10e2; 10e2 10e2 10e2])

## MooneyRivlin model
model = MooneyRivlin
u₀ = ComponentVector(C10=20e3, C01=-20e3)

## NeoHookean model
model = NeoHookean
u₀ = ComponentVector(μ=20e3)

## Isihara Model
model = Isihara
u₀ = ComponentVector(C10=20e3, C20=0.2, C01=10.0)

## Biderman Model
model = Biderman
u₀ = ComponentVector(C10=10e3, C01=100.0, C20=100.0, C30=100.0)

## JamesGreenSimpson Model
model = JamesGreenSimpson
u₀ = ComponentVector(C10=10e3, C01=1e2, C11=1e2, C20=1e2, C30=1e2)

## HainesWilson Model
model = HainesWilson
u₀ = ComponentVector(C10=1e2, C01=1e1, C11=1e1, C02=1e1, C20=1e1, C30=1e1)

## Lion Model
model = Lion
u₀ = ComponentVector(C10=1e3, C01=1e2, C50=1e3)

## HauptSedlan Model
model = HauptSedlan
u₀ = ComponentVector(C10=1e3, C01=1e1, C11=1e1, C02=1e1, C30=1e1)

## HartmannNeff Model
model = HartmannNeff
u₀ = ComponentVector(α=1e1, Ci0=[1e2], C0j=[1e1])

## Carrol Model
model = Carroll
u₀ = ComponentVector(A=1e3, B=1e1, C=1e2)

## Nunes Model
model = Nunes
u₀ = ComponentVector(C1=1e4, C2=1e3)
lb = ComponentVector(C1=0.0, C2=0.0)
ub = ComponentVector(C1=Inf, C2=Inf)

## BahremanDarijani Model
model = BahremanDarijani
u₀ = ComponentVector(A2=1e2, B2=1e2, A4=1e2, A6=1e2)

## Zhao Model
model = Zhao
u₀ = ComponentVector(C₋₁¹=1e3, C₁¹=1e3, C₂¹=1e3, C₂²=1e3)

## Knowles
model = Knowles
u₀ = ComponentVector(μ=30e3, b=6e-3, n=10000.0)
lb = ComponentVector(μ=0.0, b=0.0, n=0.0)
ub = ComponentVector(μ=Inf, b=Inf, n=Inf)

## Swanson
model = Swanson
u₀ = ComponentVector(A=[1e3, 1e3], α=[0.4, 0.3], B=[1e3, 1e2], β=[0.1, 0.1])

## Yamashita-Kawabata
model = YamashitaKawabata
u₀ = ComponentVector(C1=1e3, C2=1e3, C3=1e3, N=3.0)

## Davis-De-Thomas
model = DavisDeThomas
u₀ = ComponentVector(A=20e3, n=-0.04, C=3.5e3, k=57.6)

## Gregory - maybe?
model = Gregory
u₀ = ComponentVector(A=0.0, B=10e6, C=0.75, m=6.0, n=50.0)

## Modified Gregory - not working?
model = ModifiedGregory
u₀ = ComponentVector(A=10e3, α=10.0, M=10.0, B=20e3, β=10.0, N=10.0)

## Beda
model = Beda
u₀ = ComponentVector(C1=15e3, C2=20e3, C3=20e3, K1=30e3, α=1.0, β=0.5, ζ=1.5)
lb = ComponentVector(C1=-Inf, C2=-Inf, C3=-Inf, K1=-Inf, α=0.0, β=0.0, ζ=1.0)
ub = ComponentVector(C1=Inf, C2=Inf, C3=Inf, K1=Inf, α=Inf, β=1.0, ζ=Inf)

## Amin
model = Amin
u₀ = ComponentVector(C1=10e3, C2=4.5, C3=4.5, C4=4.5, N=1.5, M=1.5)

## LopezPamies
model = LopezPamies
u₀ = ComponentVector(α=[10.0, 1.00], μ=[10e3, 10e3])

## GenYeoh
model = GenYeoh
u₀ = ComponentVector(K1=10e3, K2=5e3, K3=20e3, m=1.0, p=1.0, q=1.0)

## VerondaWestmann
model = VerondaWestmann
u₀ = ComponentVector(C1=20e4, C2=20.0, α=-0.0)

## FungDemiray
model = FungDemiray
u₀ = ComponentVector(μ=30e3, b=-0.005)

## Vito
model = Vito
u₀ = ComponentVector(μ=20e3, b=0.005, α=0.4)

## HumphreyYin
model = HumphreyYin
u₀ = ComponentVector(C1=1e6, C2=1e-2)

## ModifiedYeoh
model = ModifiedYeoh
u₀ = ComponentVector(C10=20e3, C20=-3e1, C30=1.0, α=100.0, β=50.0)

## MansouriDarijani
model = MansouriDarijani
u₀ = ComponentVector(A1=240e3, m1=0.024, B1=0.800e6, n1=0.049)

## GentThomas
model = GentThomas
u₀ = ComponentVector(C1=20e3, C2=30e5)

## HossMarczakI
model = HossMarczakI
u₀ = ComponentVector(α=10e3, β=10e3, μ=10e3, b=1.0, n=1.0)
lb = ComponentVector(α=-Inf, β=-Inf, μ=-Inf, b=0.0, n=0.0)
ub = ComponentVector(α=Inf, β=Inf, μ=Inf, b=Inf, n=Inf)

## HossMarczakII
model = HossMarczakII
u₀ = ComponentVector(α=10e3, β=1.0, μ=10e4, b=1.0, n=0.5, C2=60e3)
lb = ComponentVector(α=-Inf, β=-Inf, μ=-Inf, b=0.0, n=0.0, C2=-Inf)
ub = ComponentVector(α=Inf, β=Inf, μ=Inf, b=Inf, n=Inf, C2=Inf)

## ExpLn
model = ExpLn
u₀ = ComponentVector(A=10e3, a=0.5, b=0.5)

## Gent
model = Gent
u₀ = ComponentVector(μ=17e3, Jₘ=1.0)
Jₘ_min = sum(collect.(data.λ⃗)[end] .^ (2)) - 3
lb = ComponentVector(μ=1.0, Jₘ=Jₘ_min)
ub = ComponentVector(μ=Inf, Jₘ=Inf)

## Killian
model = Killian
u₀ = ComponentVector(μ=20e3,)

## ArrudaBoyce
model = ArrudaBoyce
u₀ = ComponentArray(μ=243221.2739932059, N=18.29074828993436)

##  Generic Model Parameter Identification
HEProblem = HyperelasticProblem(data, model, u₀, [])
HEProblem = HyperelasticProblem(data, model, u₀, [], lb=lb, ub=ub)
sol = solve(HEProblem, BFGS())

## Sussman-Bathe Model
# W = SussmanBathe
# p = ComponentVector(s⃗ = data.s⃗, λ⃗ = data.λ⃗, k = 4)

##
s⃗ = s⃗̂(model, sol.u, collect.(data.λ⃗))
s₁ = getindex.(s⃗, 1)
# plot(getindex.(data.λ⃗, 1), (s₁ .- getindex.(data.s⃗, 1)) ./ getindex.(data.s⃗, 1) .* 100, label="Error [%]")
p = plot(getindex.(data.λ⃗, 1), s₁, label="Predicted")
plot!(p, getindex.(data.λ⃗, 1), getindex.(data.s⃗, 1), label="True")

##
# e = range(-2.0, 2.0, length=201)
# λ = exp.(e)
# τ(λ) = exp(2log(λ)) - exp(-log(λ)) + 0.5 * (1 - cos(2π * log(λ))) * (0 ≤ log(λ) ≤ 1)
# data = uniaxial_data(τ.(λ), λ)
# W = SussmanBathe(getindex.(data.s⃗, 1), getindex.(data.λ⃗, 1), 5)

# W([2.0, 1 / sqrt(2), 1 / sqrt(2)])
# ## Make Predictions with Fitted Parameters
# ∂W(x⃗) = AD.gradient(AD.FiniteDifferencesBackend(), x -> sum(W(x)), x⃗)[1]
# ∂W([2.0, 1 / sqrt(2), 1 / sqrt(2)])
# s = ∂W.(collect.(data.λ⃗))

# s₁ = getindex.(s, 1) - getindex.(s, 3)
# p = plot(getindex.(data.λ⃗, 1), (s₁ .- getindex.(data.s⃗, 1)) ./ getindex.(data.s⃗, 1) .* 100, label="Error [%]")

# p = plot(getindex.(data.λ⃗, 1), s₁, label="Predicted")
# plot!(p, getindex.(data.λ⃗, 1), getindex.(data.s⃗, 1), label="True")