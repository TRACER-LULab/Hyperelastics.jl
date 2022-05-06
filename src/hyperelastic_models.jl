module HyperelasticModels

using Tullio
using SpecialFunctions
export GeneralizedMooneyRivlin, GeneralizedDarijaniNaghdabadi, GeneralizedBeda, MooneyRivlin, NeoHookean, Gent, Biderman, Isihara, JamesGreenSimpson, Lion, Yeoh, HauptSedlan, HartmannNeff, HainesWilson, Carroll, BahremanDarijani, Zhao, Knowles, Swanson, YamashitaKawabata, DavisDeThomas, Gregory, ModifiedGregory, Beda, Amin, LopezPamies, GenYeoh, HartSmith, VerondaWestmann, FungDemiray, Vito, ModifiedYeoh, Martins, ChevalierMarco, GornetDesmorat, MansouriDarijani, GentThomas, Alexander, LambertDianiRey, HossMarczakI, HossMarczakII, ExpLn, Kilian, VanDeWaals, TakamizawaHayashi, YeohFleming, PucciSaccomandi, HorganSaccomandi, Beatty, HorganMurphy, ArrudaBoyce, Ogden, EdwardVilgis

export ValanisLandel, PengLandel, Ogden, Attard, Shariff, ArmanNarooei

"""
General Mooney Rivlin

Parameters: [C]

Model: ``∑\\limits_{i,j = 0}^{N,M} Cᵢⱼ(I₁-3)ⁱ(I₂-3)ʲ``
"""
function GeneralizedMooneyRivlin((; C))
    function (λ⃗)
        I1 = I₁(λ⃗)
        I2 = I₂(λ⃗)
        @tullio W := C[j, i] * (I1 - 3)^(i - 1) * (I2 - 3)^(j - 1)
        return W
    end
end

function GeneralizedDarijaniNaghdabadi((; A, B, m, n))
    (λ⃗) -> sum(A .* (λ⃗ .^ m .- 1) + B .* (λ⃗ .^ (-n) .- 1))
end

function GeneralizedBeda((; C, K, α, β))
    function (λ⃗)
        # @tullio W:= C[i]/α[i]*(I₁(λ⃗)-3)^α[i]+K[j]/β[j]*(I₂(λ⃗)-3)^β[j]
        W1 = C ./ α .* (I₁(λ⃗) - 3) .^ α |> sum
        W2 = K ./ β .* (I₂(λ⃗) - 3) .^ β |> sum
        return W1 + W2
    end
end

function MooneyRivlin((; C10, C01))
    W = GeneralizedMooneyRivlin(ComponentVector(
        C=[
            0.0 C10
            C01 0.0
        ]))
    (λ⃗) -> W(λ⃗)
end

function NeoHookean((; μ))
    W = GeneralizedMooneyRivlin(ComponentVector(
        C=[
        0.0 μ
    ]))
    (λ⃗) -> W(λ⃗)
end

function Isihara((; C10, C20, C01))
    W = GeneralizedMooneyRivlin(ComponentVector(
        C=[
            0.0 C10 C20
            C01 0.0 0.0
        ]))
    (λ⃗) -> W(λ⃗)
end

function Biderman((; C10, C01, C20, C30))
    W = GeneralizedMooneyRivlin(ComponentVector(
        C=[
            0.0 C10 C20 C30
            C01 0.0 0.0 0.0
        ]
    ))
    (λ⃗) -> W(λ⃗)
end

function JamesGreenSimpson((; C10, C01, C11, C20, C30))
    W = GeneralizedMooneyRivlin(ComponentVector(
        C=[
            0.0 C10 C20 C30
            C01 C11 0.0 0.0
        ]))
    (λ⃗) -> W(λ⃗)
end

function HainesWilson((; C10, C01, C11, C02, C20, C30))
    W = GeneralizedMooneyRivlin(ComponentVector(
        C=[
            0.0 C10 C20 C30
            C01 C11 0.0 0.0
            C02 0.0 0.0 0.0
        ]))
    (λ⃗) -> W(λ⃗)
end

function Yeoh((; C10, C20, C30))
    W = GeneralizedMooneyRivlin(ComponentVector(
        C=[
        0.0 C10 C20 C30
    ]))
    (λ⃗) -> W(λ⃗)
end

function Lion((; C10, C01, C50))
    W = GeneralizedMooneyRivlin(ComponentVector(
        C=[
            0.0 C10 0.0 0.0 0.0 C50
            C01 0.0 0.0 0.0 0.0 0.0
        ]))
    (λ⃗) -> W(λ⃗)
end

function HauptSedlan((; C10, C01, C11, C02, C30))
    W = GeneralizedMooneyRivlin(ComponentVector(
        C=[
            0.0 C10 0.0 C30
            C01 C11 0.0 0.0
            C02 0.0 0.0 0.0
        ]))
    (λ⃗) -> W(λ⃗)
end

function HartmannNeff((; α, Ci0, C0j))
    function f(λ⃗)
        @tullio ∑ = Ci0[i] * (I₁(λ⃗) - 3)^i + C0j[j] * (I₂(λ⃗)^(3 / 2) - 3sqrt(3))^j
        α * (I₁(λ⃗)^3 - 3) + ∑
    end
end

function Carroll((; A, B, C))
    (λ⃗) -> A * I₁(λ⃗) + B * I₁(λ⃗)^4 + C * I₂(λ⃗)^(1 / 2)
end

## Only developed for simple shear deformation.
# function Nunes((; C1, C2))
# (λ⃗) -> C1 * (I₁(λ⃗) - 3) + 4 / 3 * C2 * (I₂(λ⃗) - 3)^(3 / 4)
# (λ⃗) -> 1/2*(C1*(I₁(λ⃗)-3)+C2*(I₂(λ⃗)-3)^(3/4))
# end

function BahremanDarijani((; A2, B2, A4, A6))
    W = GeneralizedDarijaniNaghdabadi(
        ComponentVector(
            A=[0.0, A2, 0.0, A4, 0.0, A6],
            B=[0.0, B2],
            m=[0.0, 2.0, 0.0, 4.0, 0.0, 6.0],
            n=[0.0, 2.0])
    )
    (λ⃗) -> W(λ⃗)
end

function Zhao((; C₋₁¹, C₁¹, C₂¹, C₂²))
    (λ⃗) -> C₋₁¹ * (I₂(λ⃗) - 3) + C₁¹ * (I₁(λ⃗) - 3) + C₂¹ * (I₁(λ⃗)^2 - 2I₂(λ⃗) - 3) + C₂² * (I₁(λ⃗)^2 - 2I₂(λ⃗) - 3)^2
end

## Table 2
function Knowles((; μ, b, n))
    (λ⃗) -> μ / (2b) * ((1 + (b / n) * (I₁(λ⃗) - 3))^n - 1)
end

# Article Requested
function Swanson((; A, α, B, β))
    (λ⃗) -> @tullio _ := 3 / 2 * (A[i] / (1 + α[i]) * (I₁(λ⃗) / 3)^(1 + α[i]) + B[i] / (1 + β[i]) * (I₂(λ⃗) / 3)^(1 + β[i]))
end

# Original article in Japanese
function YamashitaKawabata((; C1, C2, C3, N))
    (λ⃗) -> C1 * (I₁(λ⃗) - 3) + C2 * (I₂(λ⃗) - 3) + C3 / (N + 1) * (I₁(λ⃗) - 3)^(N + 1)
end

# Article Requested
function DavisDeThomas((; A, n, C, k))
    (λ⃗) -> A / (2 * (1 - n / 2)) * (I₁(λ⃗) - 3 + C^2)^(1 - n / 2) + k * (I₁(λ⃗) - 3)^2
end

# Article Requested
function Gregory((; A, B, C, m, n))
    (λ⃗) -> A / (2 - n) * (I₁(λ⃗) - 3 + C^2)^(1 - n / 2) + B / (2 + m) * (I₁(λ⃗) - 3 + C^2)^(1 + m / 2)
end

# Proposed in 85 Model review
function ModifiedGregory((; A, α, M, B, β, N))
    (λ⃗) -> A / (1 + α) * (I₁(λ⃗) - 3 + M^2)^(1 + α) + B / (1 + β) * (I₁(λ⃗) - 3 + N^2)^(1 + β)
end

# Added general form of the Beda model
function Beda((; C1, C2, C3, K1, α, β, ζ))
    W = GeneralizedBeda(ComponentVector(
        C=[C1, C2, C3],
        K=[K1],
        α=[α, 1.0, ζ],
        β=[β]
    )
    )
    (λ⃗) -> W(λ⃗)
    # (λ⃗) -> C1 / α * (I₁(λ⃗) - 3)^(α) + C2 * (I₁(λ⃗) - 3) + C3 / ζ * (I₁(λ⃗) - 3)^(ζ) + K1 / β * (I₂(λ⃗) - 3)^β
end

function Amin((; C1, C2, C3, C4, N, M))
    (λ⃗) -> C1 * (I₁(λ⃗) - 3) + C2 / (N + 1) * (I₁(λ⃗) - 3)^(N + 1) + C3 / (M + 1) * (I₁(λ⃗) - 3)^(M + 1) + C4 * (I₂(λ⃗) - 3)
end

# Article Received - General form presented.
function LopezPamies((; α, μ))
    (λ⃗) -> @tullio _ := (3.0^(1 - α[i])) / (2α[i]) * μ[i] * (I₁(λ⃗)^(α[i]) - 3^(α[i]))
end

# ✔
function GenYeoh((; K1, K2, K3, m, p, q))
    (λ⃗) -> K1 * (I₁(λ⃗) - 3)^m + K2 * (I₁(λ⃗) - 3)^p + K3 * (I₁(λ⃗) - 3)^q
end
# ✔
function VerondaWestmann((; C1, C2, α))
    (λ⃗) -> C1 * (exp(α * (I₁(λ⃗) - 3)) - 1) + C2 * (I₂(λ⃗) - 3)
end
# ✔
function FungDemiray((; μ, b))
    (λ⃗) -> μ / (2 * b) * (exp(b * (I₁(λ⃗) - 3)) - 1)
end
# ✔
function Vito((; α, β, γ))
    (λ⃗) -> α * (exp(β * (I₁(λ⃗) - 3) + γ * (I₂(λ⃗) - 3)) - 1)
end

# Only applicable for fiber composites
# function HumphreyYin((; c, b, A, a))
#     (λ⃗, α) -> c * (exp(b * (I₁(λ⃗) - 3)) - 1) + A*(exp(exp(a*(α-1)^2))-1)
# end

# Requested in ILL
function ModifiedYeoh((; C10, C20, C30, α, β))
    (λ⃗) -> C10 * (I₁(λ⃗) - 3) + C20 * (I₁(λ⃗) - 3)^2 + C30 * (I₁(λ⃗) - 3)^3 + α / β * (1 - exp(-β * (I₁(λ⃗) - 3)))
end
# ✔
function MansouriDarijani((; A1, m1, B1, n1))
    (λ⃗) -> A1 * (exp(m1 * (I₁(λ⃗) - 3)) - 1) + B1 * (exp(n1 * (I₂(λ⃗) - 3)) - 1)
end
# ✔
function GentThomas((; C1, C2))
    (λ⃗) -> C1 * (I₁(λ⃗) - 3) + C2 * log(I₂(λ⃗) / 3)
end
# Suitable for low strains
function HossMarczakI((; α, β, μ, b, n))
    (λ⃗) -> α / β * (1 - exp(-β * (I₁(λ⃗) - 3))) + μ / (2b) * ((1 + b / n * (I₁(λ⃗) - 3))^n - 1)
end
# Suitable for high strains
function HossMarczakII((; α, β, μ, b, n, C2))
    (λ⃗) -> α / β * (1 - exp(-β * (I₁(λ⃗) - 3))) + μ / (2b) * ((1 + b / n * (I₁(λ⃗) - 3))^n - 1) + C2 * log(I₂(λ⃗) / 3)
end

function ExpLn((; A, a, b))
    (λ⃗) -> A * (1 / a * exp(a * (I₁(λ⃗) - 3)) + b * (I₁(λ⃗) - 2) * (1 - log(I₁(λ⃗) - 2)) - 1 / a - b)
end

#####################
###### TABLE 3 ######
#####################
# model not mentioned in original article ?????? 
# function Warner((; μ, Iₘ))
#     (λ⃗) -> -1/2*μ*Iₘlog(1-(I₁(λ⃗)-3)/(Iₘ-3))
# end 

# does not easily match the form in the paper
# function Killian((; μ, JL))
#     (λ⃗) -> -μ * JL * (log(1 - sqrt((I₁(λ⃗) - 3) / JL)) + sqrt((I₁(λ⃗) - 3) / JL))
# end
# ARticles requested -> Checked against review article from Marckmann and Verron
function VanDerWaals((; μ, λm, β, α))
    (λ⃗) -> μ * (-(λm^2 - 3) * (log(1 - θ) + θ) - 2 / 3 * α * ((I₁(λ⃗) - 3) / 2)^(3 / 2))
end

function Gent((; μ, Jₘ))
    (λ⃗) -> -μ * Jₘ / 2 * log(1 - (I₁(λ⃗) - 3) / Jₘ)
end

# With the assumption of isotropicity -> Verified with A description of arterial wall mechanics using limiting chain extensibility constitutitive models by Horgan and Saccomandi
function TakamizawaHayashi((; c, Jm))
    (λ⃗) -> -c * log(1 - ((I₁(λ⃗) - 3) / Jm)^2)
end

function YeohFleming((; A, B, C10, Im))
    (λ⃗) -> A / B * (1 - exp(-B * (I₁(λ⃗) - 3))) - C10 * (Im - 3) * log(1 - ((I₁(λ⃗) - 3) / (Im - 3)))
end

## Not a real model => Not referenced in the Gent paper cited
# function Gent3Parameters((;μ, Jₘ, α))
#     (λ⃗) -> 
# end

function PucciSaccomandi((; K, μ, Jₘ))
    (λ⃗) -> K * log(I₂(λ⃗) / 3) - μ * Jₘ / 2 * log(1 - (I₁(λ⃗) - 3) / Jₘ)
end

# Originally from CONSTITUTIVE MODELS FOR ATACTIC ELASTOMERS
"""
Horgan Saccomandi Model

Parameters: μ, J

Model: ``-\\frac{μJ}{2}\\log\\bigg(\\frac{J^3-J^2I₁+JI₂-1}{(J-1)^3}\\bigg)```
"""
function HorganSaccomandi((; μ, J))
    (λ⃗) -> -μ * J / 2 * log((J^3 - J^2 * I₁(λ⃗) + J * I₂(λ⃗) - 1) / (J - 1)^3)
end

"""
Beatty Model

Parameters: G₀, Iₘ

Model: ``-\\frac{G₀Iₘ(Iₘ-3)}{2(2Iₘ-3)}\\log\\bigg(\\frac{1-\\frac{I₁-3}{Iₘ-3}}{1+\\frac{I₁-3}{Iₘ}} \\bigg)``
"""
function Beatty((; G₀, Iₘ))
    (λ⃗) -> -G₀ * Iₘ * (Iₘ - 3) / 2 / (2Iₘ - 3) * log((1 - (I₁(λ⃗) - 3) / (Iₘ - 3)) / (1 + (I₁(λ⃗) - 3) / (Iₘ)))
end

"""
Horgan Murphy Model

Parameters: μ, Jₘ, c

Model: ``-\\frac{2μJₘ}{c^2}\\log\\bigg(1-\\frac{λ₁ᶜ+λ₂ᶜ+λ₃ᶜ-3}{Jₘ})``
"""
function HorganMurphy((; μ, Jₘ, c))
    (λ⃗) -> -2 * μ * Jₘ / c^2 * log(1 - (sum(λ⃗ .^ c) - 3) / Jₘ)
end

########################
###########33 TABLE 4
########################
"""
Valanis-Landel

Parameters: μ

Model: ``2μ∑₁³(λᵢ(\\log\\lambda_i -1))
"""
function ValanisLandel((; μ))
    (λ⃗) -> 2 * μ * sum(λ⃗ * (log.(λ⃗) - 1))
end

"""
Peng - Landel

Parameters: E

Model: ``E∑₁³\\bigg[λᵢ - 1 - \\log(λᵢ) - \\frac{1}{6}\\log(λᵢ)² + \\frac{1}{18}\\log(λᵢ)³-\\frac{1}{216}\\log(λᵢ)⁴]
"""
function PengLandel((; E))
    (λ⃗) -> sum(@. λ⃗ - 1 - log(λ⃗) - 1 / 6 * log(λ⃗)^2 + 1 / 18 * log(λ⃗)^3 - 1 / 216 * log(λ⃗)^4) * E
end

"""
Ogden

Parameters: μ⃗, α⃗

Model: ``∑₁ᴺ \\frac{μᵢ}{αᵢ}(λ₁^αᵢ+λ₂^αᵢ+λ₃^αᵢ-3)``
"""
function Ogden((; μ, α))
    (λ⃗) -> @tullio _ := μ[i] / α[i] * (sum(λ⃗ .^ α[i]) - 3)
end

"""
Attard

Parameters: A⃗, B⃗

Model: ``∑₁ᴺ\\frac{Aᵢ}{2i}(λ₁^{2i}+λ₂^{2i}+λ₃^{2i}-3) + \\frac{Bᵢ}{2i}(λ₁^{-2i}+λ₂^{-2i}+λ₃^{-2i}-3)``
"""
function Attard((; A, B))
    (λ⃗) -> @tullio _ := A[i] / 2 / i * (sum(λ⃗ .^ (2i)) - 3) + B[i] / 2 / i * (sum(λ⃗ .^ (-2i)) - 3)
end

"""
Shariff

Parameters: E, α₁, α₂, α₃, α₄, α₅

Model: ``E*∑ᵢ∑ⱼαⱼΦⱼ(λᵢ)``
"""
function Shariff((; E, α))
    ϕ = []
    c(j, r) = factorial(j) / factorial(r) / factorial(j - r)
    for j in eachindex(α)
        if j == 0
            push!(ϕ, x -> ln(x)^2 / 3)
        elseif j == 1
            push!(ϕ, x -> -exp(1) * expinti(-1) + exp(1) * expinti(-x) + x - 2log(x) - 1)
        elseif j == 2
            push!(ϕ, x -> (expinti(x) - expinti(1)) / exp(1) - x + 1)
        elseif j == 3
            push!(ϕ, x -> -1 / (0.6 * x^(0.6)) + 3 / (1.6 * x^(1.6)) - 3 / (2.6 * x^(2.6)) + 1 / (5.6 * x^(5.6)) + 107200 / 139776)
        else
            push!(ϕ, x -> (-1)^(j - 1) * log(x) + (-1)^(j - 1) * sum(r -> (-1)^r * c(j - 1, r) * x^r / r, range(1, j - 1)) - (-1)^(j - 1) * sum(r -> (-1)^r * c(j - 1, r) / r, range(1, j - 1)))
        end
    end
    (λ⃗) -> E * (@tullio _ := ϕ[i](λ⃗[j]))
end

# Article requested
"""
Armna - Narooei

Parameters: A⃗, B⃗, m⃗, n⃗, α⃗, β⃗

Model: ``∑ᵢᴺ Aᵢ[exp(mᵢ(λ₁^{αᵢ}+λ₂^{αᵢ}+λ₃^{αᵢ}-3))-1] + Bᵢ[exp(nᵢ(λ₁^{-βᵢ}+λ₂^{-βᵢ}+λ₃^{-βᵢ}-3))-1] ``
"""
function ArmanNarooei((; A, B, m, n, α, β))
    (λ⃗) -> @tullio _ := A[i] * (exp(m[i] * (sum(λ⃗ .^ α[i]) - 3)) - 1) + B[i] * (exp(n[i] * (sum(λ⃗ .^ (-β[i])) - 3)) - 1)
end
# ---------------------------------------------------- #


function ArrudaBoyce((; μ, N))
    (λ⃗) -> μ * (0.5 * (I₁(λ⃗) - 3) + 1 / 20 / N * (I₁(λ⃗)^2 - 9) + 11 / 1050 / N^2 * (I₁(λ⃗) - 27) + 19 / 7000 / N^3 * (I₁(λ⃗)^4 - 81) + 519 / 673750 / N^4 * (I₁(λ⃗)^5 - 243))
end



function EdwardVilgis((; Ns, Nc, α, η))
    function W(λ⃗)
        A = 0.5 * Nc * ((1 - α^2) * I₁(λ⃗) / (1 - α^2 * I₁(λ⃗)) + log(1 - α^2 * I₁(λ⃗)))
        B = sum(i -> (1 + η) * (1 - α^2) * λ⃗[i] / (1 - η * λ⃗[i]^2) / (1 - α^2 * I₁(λ⃗)) + log(1 + η * λ⃗[i]^2), 1:3)
        B = 0.5 * Ns * (B + log(1 - α^2 * I₁(λ⃗)))
        W = A + B
        return W
    end
end
end