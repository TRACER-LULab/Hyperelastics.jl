module HyperelasticModels

using Tullio
using SpecialFunctions
using ComponentArrays

export GeneralMooneyRivlin, GeneralDarijaniNaghdabadi, GeneralBeda, MooneyRivlin, NeoHookean, Gent, Biderman, Isihara, JamesGreenSimpson, Lion, Yeoh, HauptSedlan, HartmannNeff, HainesWilson, Carroll, BahremanDarijani, Zhao, Knowles, Swanson, YamashitaKawabata, DavisDeThomas, Gregory, ModifiedGregory, Beda, Amin, LopezPamies, GenYeoh, HartSmith, VerondaWestmann, FungDemiray, Vito, ModifiedYeoh, Martins, ChevalierMarco, GornetDesmorat, MansouriDarijani, GentThomas, Alexander, LambertDianiRey, HossMarczakI, HossMarczakII, ExpLn, Kilian, VanDerWaals, TakamizawaHayashi, YeohFleming, PucciSaccomandi, HorganSaccomandi, Beatty, HorganMurphy, ArrudaBoyce, Ogden, EdwardVilgis

export ValanisLandel, PengLandel, Ogden, Attard, Shariff, ArmanNarooei

include("BasicDefinitions.jl")
# """
# First stretch invariant - Currently requires the addition of 5 times the machine precision to allow AD to work correctly

# ``I_1 = \\lambda_1^2+\\lambda_2^2+\\lambda_3^2 + 5\\varepsilon``
# """
# function I₁(λ⃗)
#     sum(λ⃗ .^ 2) + 5eps(Float64)
# end

# """
# Second Stretch invariant

# ``I_2 = \\lambda_1^{-2}+\\lambda_2^{-2}+\\lambda_3^{-2}``
# """
# function I₂(λ⃗)
#     sum(λ⃗ .^ (-2)) + 5eps(Float64)
# end

# """
# Third Stretch invariant

# ``I_3 = (\\lambda_1\\lambda_\\lamdba_3)^2``
# """
# function I₃(λ⃗)
#     prod(λ⃗)^2
# end

# """
# Volumetric Stretch

# ``J = \\lambda_1\\lambda_2\\lambda_3``
# """
# function J(λ⃗)
#     prod(λ⃗)
# end
"""
General Mooney Rivlin

Parameters: [C]

Model: ``\\sum\\limits_{i,j = 0}^{N,M} C_{i,j}(I_1-3)^i(I_2-3)^j``
"""
function GeneralMooneyRivlin((; C))
    function (λ⃗)
        I1 = I₁(λ⃗)
        I2 = I₂(λ⃗)
        @tullio W := C[j, i] * (I1 - 3)^(i - 1) * (I2 - 3)^(j - 1)
        return W
    end
end

"""
General Darijani Naghdabadi

Parameters: A⃗, B⃗, m⃗, n⃗

Model: ``\\sum\\limits_{i = 1}{3}\\sum\\limits_{j=0}^{N} A_j (\\lambda_i^{m_j}-1) + B_j(\\lambda_i^{-n_j}-1)``
"""
function GeneralDarijaniNaghdabadi((; A⃗, B⃗, m⃗, n⃗))
    (λ⃗) -> sum(A⃗ .* (λ⃗ .^ m⃗ .- 1) + B⃗ .* (λ⃗ .^ (-n⃗) .- 1))
end

"""
General Beda

Parameters: C, K, α, β

Model: ``\\sum\\limits_{i = 1}^{N}\\frac{C_i}{\\alpha_i}(I_1-3)^{\\alpha_i} + \\sum\\limits_{j=1}^{M}\\frac{K_j}{\\beta_j}(I_2-3)^{\\beta_j}``
"""
function GeneralBeda((; C, K, α, β))
    function (λ⃗)
        # @tullio W:= C[i]/α[i]*(I₁(λ⃗)-3)^α[i]+K[j]/β[j]*(I₂(λ⃗)-3)^β[j]
        W1 = C ./ α .* (I₁(λ⃗) - 3) .^ α |> sum
        W2 = K ./ β .* (I₂(λ⃗) - 3) .^ β |> sum
        return W1 + W2
    end
end

"""
Mooney Rivlin Model

Parameters: C01, C10

Model: ``C_{10}(I_1-3)+C_{01}(I_2-3)``
"""
function MooneyRivlin((; C10, C01))
    W = GeneralMooneyRivlin(ComponentVector(
        C=[
            0.0 C10
            C01 0.0
        ]))
    (λ⃗) -> W(λ⃗)
end

"""
NeoHookean

Parameters: μ

Model: ``\\frac{\\mu}{2}(I_1-3)``
"""
function NeoHookean((; μ))
    W = GeneralMooneyRivlin(ComponentVector(
        C=[
        0.0 μ / 2
    ]))
    (λ⃗) -> W(λ⃗)
end

"""
Isihara 

Parameters: C10, C20, C01

Model: ``\\sum\\limits_{i,j=0}^{2, 1}C_{i,j}(I_1-3)^i(I_2-3)^j``
"""
function Isihara((; C10, C20, C01))
    W = GeneralMooneyRivlin(ComponentVector(
        C=[
            0.0 C10 C20
            C01 0.0 0.0
        ]))
    (λ⃗) -> W(λ⃗)
end

"""
Biderman 

Parameters: C10, C01, C20, C30

Model: ``\\sum\\limits_{i,j=0}^{3, 1}C_{i,j}(I_1-3)^i(I_2-3)^j``
"""
function Biderman((; C10, C01, C20, C30))
    W = GeneralMooneyRivlin(ComponentVector(
        C=[
            0.0 C10 C20 C30
            C01 0.0 0.0 0.0
        ]
    ))
    (λ⃗) -> W(λ⃗)
end

"""
James-Green-Simpson 

Parameters: C10, C01, C11, C20, C30

Model: ``\\sum\\limits_{i,j=0}^{3, 1}C_{i,j}(I_1-3)^i(I_2-3)^j``
"""
function JamesGreenSimpson((; C10, C01, C11, C20, C30))
    W = GeneralMooneyRivlin(ComponentVector(
        C=[
            0.0 C10 C20 C30
            C01 C11 0.0 0.0
        ]))
    (λ⃗) -> W(λ⃗)
end

"""
Haines-Wilson

Parameters: C10, C01, C11, C02, C20, C30

Model: ``\\sum\\limits_{i,j=0}^{3, 2}C_{i,j}(I_1-3)^i(I_2-3)^j``
"""
function HainesWilson((; C10, C01, C11, C02, C20, C30))
    W = GeneralMooneyRivlin(ComponentVector(
        C=[
            0.0 C10 C20 C30
            C01 C11 0.0 0.0
            C02 0.0 0.0 0.0
        ]))
    (λ⃗) -> W(λ⃗)
end

"""
Yeoh

Parameters: C10, C20, C30

Model: ``\\sum\\limits_{i,j=0}^{3, 0}C_{i,j}(I_1-3)^i(I_2-3)^j``
"""
function Yeoh((; C10, C20, C30))
    W = GeneralMooneyRivlin(ComponentVector(
        C=[
        0.0 C10 C20 C30
    ]))
    (λ⃗) -> W(λ⃗)
end

"""
Lion

Parameters: C10, C01, C50

Model: ``\\sum\\limits_{i,j=0}^{5,1}C_{i,j}(I_1-3)^i(I_2-3)^j``
"""
function Lion((; C10, C01, C50))
    W = GeneralMooneyRivlin(ComponentVector(
        C=[
            0.0 C10 0.0 0.0 0.0 C50
            C01 0.0 0.0 0.0 0.0 0.0
        ]))
    (λ⃗) -> W(λ⃗)
end

"""
Haupt Sedlan

Parameters: C10, C01, C11, C02, C30

Model: ``\\sum\\limits_{i,j=0}^{3, 2}C_{i,j}(I_1-3)^i(I_2-3)^j``
"""
function HauptSedlan((; C10, C01, C11, C02, C30))
    W = GeneralMooneyRivlin(ComponentVector(
        C=[
            0.0 C10 0.0 C30
            C01 C11 0.0 0.0
            C02 0.0 0.0 0.0
        ]))
    (λ⃗) -> W(λ⃗)
end

"""
Hartmann-Neff

Parameters: α, Ci0, C0j

Model: ``\\sum\\limits_{i,j=0}^{M,N}C_{i,0}(I_1-3)^i -3\\sqrt{3}^j+\\alpha(I_1-3)``
"""
function HartmannNeff((; α, Ci0, C0j))
    function f(λ⃗)
        @tullio ∑ = Ci0[i] * (I₁(λ⃗) - 3)^i + C0j[j] * (I₂(λ⃗)^(3 / 2) - 3sqrt(3))^j
        α * (I₁(λ⃗)^3 - 3) + ∑
    end
end

"""
Carroll

Parameters: A, B, C 

Model: ``AI_1+BI_1^4+C\\sqrt{I_2}``
"""
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
"""
Zhao

Parameters: C₋₁¹,, C₁¹, C₂¹, C₂²

Model: ``C_{-1}^1*(I_2-3)+C_{1}^{1}(I_1-3)+C_{2}^{1}(I_1^2-2I_2-3)+C_{2}^{2}(I_1^2-2I_2-3)^2``
"""
function Zhao((; C₋₁¹, C₁¹, C₂¹, C₂²))
    (λ⃗) -> C₋₁¹ * (I₂(λ⃗) - 3) + C₁¹ * (I₁(λ⃗) - 3) + C₂¹ * (I₁(λ⃗)^2 - 2I₂(λ⃗) - 3) + C₂² * (I₁(λ⃗)^2 - 2I₂(λ⃗) - 3)^2
end

## Table 2
"""
Knowles

Parameters: μ, b, n

Model: ``\\frac{\\mu}{2b}((1+\\frac{b}{n}(I_1-3))^n-1)``
"""
function Knowles((; μ, b, n))
    (λ⃗) -> μ / (2b) * ((1 + (b / n) * (I₁(λ⃗) - 3))^n - 1)
end

# Article Requested
"""
Swanson

Parameters: A, α, B, β

Model: ``\\sum\\limits_{i=1}^{N} \\frac{3}{2}(\\frac{A_i}{1+\\alpha_i}(\\frac{I_1}{3})^{1+\\alpha_i}+\\frac{B_i}{1+\\beta_i}(\\frac{I_2}{3})^{1+\\beta_i}``
"""
function Swanson((; A, α, B, β))
    (λ⃗) -> @tullio _ := 3 / 2 * (A[i] / (1 + α[i]) * (I₁(λ⃗) / 3)^(1 + α[i]) + B[i] / (1 + β[i]) * (I₂(λ⃗) / 3)^(1 + β[i]))
end

# Original article in Japanese
"""
Yamashita-Kawabata

Parameters: C1, C2, C3, N

Model: ``C_1(I_1-3)+C_2(I_2-3)+\\frac{C_3}{N+1}(I_1-3)^{N+1}``
"""
function YamashitaKawabata((; C1, C2, C3, N))
    (λ⃗) -> C1 * (I₁(λ⃗) - 3) + C2 * (I₂(λ⃗) - 3) + C3 / (N + 1) * (I₁(λ⃗) - 3)^(N + 1)
end

# Article Requested
"""
Davis-DeThomas

Parameters: A, n, C, k

Model: ``\\frac{A}{2(1-\\frac{n}{2})}(I_1-3+C^2)^{1-\\frac{n}{2}}+k(I_1-3)^2``
"""
function DavisDeThomas((; A, n, C, k))
    (λ⃗) -> A / (2 * (1 - n / 2)) * (I₁(λ⃗) - 3 + C^2)^(1 - n / 2) + k * (I₁(λ⃗) - 3)^2
end

# Article Requested
"""
Gregory

Parameters: A, B, C, m, n

Model: ``\\frac{A}{2-n}(I_1-3+C^2)^{1-\\frac{n}{2}}+\\frac{B}{2+m}(I_1-3+C^2)^{1+\\frac{m}{2}}``
"""
function Gregory((; A, B, C, m, n))
    (λ⃗) -> A / (2 - n) * (I₁(λ⃗) - 3 + C^2)^(1 - n / 2) + B / (2 + m) * (I₁(λ⃗) - 3 + C^2)^(1 + m / 2)
end

# Proposed in 85 Model review
"""
Modified Gregory

Parameters: A, α, M, B, β, N

Model: ``\\frac{A}{1+\\alpha}(I_1-3+M^2)^{1+\\alpha}+\\frac{B}{1+\\beta}(I_1-3+N^2)^{1+\\beta}``
"""
function ModifiedGregory((; A, α, M, B, β, N))
    (λ⃗) -> A / (1 + α) * (I₁(λ⃗) - 3 + M^2)^(1 + α) + B / (1 + β) * (I₁(λ⃗) - 3 + N^2)^(1 + β)
end

# Added general form of the Beda model
"""
Beda

Parameters: C1, C2, C3, K1, α, β, ζ

Model: ``\\frac{C_1}{\\alpha}(I_1-3)^{\\alpha}+C_2(I_1-3)+\\frac{C_3}{\\zeta}(I_1-3)^{\\zeta}+\\frac{K_1}{\\beta}(I_2-3)^\\beta``
"""
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

"""
Amin

Parameters: C1, C2, C3, C4, N, M

Model: ``C_1 (I_1 - 3) + \\frac{C_2}{N + 1} (I_1 - 3)^{N + 1} + \\frac{C_3}{M + 1} (I_1 - 3)^{M + 1} + C_4 (I_2 - 3)``
"""
function Amin((; C1, C2, C3, C4, N, M))
    (λ⃗) -> C1 * (I₁(λ⃗) - 3) + C2 / (N + 1) * (I₁(λ⃗) - 3)^(N + 1) + C3 / (M + 1) * (I₁(λ⃗) - 3)^(M + 1) + C4 * (I₂(λ⃗) - 3)
end

# Article Received - General form presented.
"""
Lopez-Pamies

Parameters: α⃗, μ⃗

Model: ``\\frac{3.0^{1 - \\alpha_i}}{2\\alpha_i} \\mu_i (I_1^{\\alpha_i} - 3^{\\alpha_i})``
"""
function LopezPamies((; α⃗, μ⃗))
    (λ⃗) -> @tullio _ := (3.0^(1 - α⃗[i])) / (2α⃗[i]) * μ⃗[i] * (I₁(λ⃗)^(α⃗[i]) - 3^(α⃗[i]))
end

# ✔
"""
GenYeoh

Parameters: K1, K2, K3, m, p, q

Model: ``K_1 (I_1 - 3)^m + K_2 * (I_1 - 3)^p + K_3 * (I_1 - 3)^q``
"""
function GenYeoh((; K1, K2, K3, m, p, q))
    (λ⃗) -> K1 * (I₁(λ⃗) - 3)^m + K2 * (I₁(λ⃗) - 3)^p + K3 * (I₁(λ⃗) - 3)^q
end
# ✔
"""
Veronda-Westmann

Parameters: C1, C2, α

Model: ``C_1 (\\exp(\\alpha(I_1 - 3)) - 1) + C_2 (I_2 - 3)``
"""
function VerondaWestmann((; C1, C2, α))
    (λ⃗) -> C1 * (exp(α * (I₁(λ⃗) - 3)) - 1) + C2 * (I₂(λ⃗) - 3)
end
# ✔
"""
Fung-Demiray

Parameters: μ, b

Model: ``\\frac{\\mu}{2 * b} (\\exp(b(I_1 - 3)) - 1)``
"""
function FungDemiray((; μ, b))
    (λ⃗) -> μ / (2 * b) * (exp(b * (I₁(λ⃗) - 3)) - 1)
end
# ✔
"""
Vito

Parameters: α, β, γ

Model: ``\\alpha (\\exp\\bigg(\\beta (I_1 - 3)\\bigg) + \\gamma * (I_2 - 3)) - 1)``
"""
function Vito((; α, β, γ))
    (λ⃗) -> α * (exp(β * (I₁(λ⃗) - 3) + γ * (I₂(λ⃗) - 3)) - 1)
end

# Only applicable for fiber composites
# function HumphreyYin((; c, b, A, a))
#     (λ⃗, α) -> c * (exp(b * (I₁(λ⃗) - 3)) - 1) + A*(exp(exp(a*(α-1)^2))-1)
# end

# Requested in ILL
"""
Modified Yeoh

Parameters: C10, C20, C30, α, β

Model: ``C_{10} * (I_1 - 3) + C_{20} * (I_1 - 3)^2 + C_{30} * (I_1 - 3)^3 + \\alpha / \\beta * (1 - \\exp{-\\beta * (I_1 - 3)})``
"""
function ModifiedYeoh((; C10, C20, C30, α, β))
    (λ⃗) -> C10 * (I₁(λ⃗) - 3) + C20 * (I₁(λ⃗) - 3)^2 + C30 * (I₁(λ⃗) - 3)^3 + α / β * (1 - exp(-β * (I₁(λ⃗) - 3)))
end
# ✔
"""
Mansouri-Darijani

Parameters: A1, m1, B1, n1

Model: ``A_1\\exp{m_1(I_1-3)-1}+B_1\\exp{n_1(I_2-3)-1}``
"""
function MansouriDarijani((; A1, m1, B1, n1))
    (λ⃗) -> A1 * (exp(m1 * (I₁(λ⃗) - 3)) - 1) + B1 * (exp(n1 * (I₂(λ⃗) - 3)) - 1)
end
# ✔
"""
Gent Thomas

Paramters: C1, C2

Model: ``C_1(I_1-3)+C_2\\log(\\frac{I_2}{3})``
"""
function GentThomas((; C1, C2))
    (λ⃗) -> C1 * (I₁(λ⃗) - 3) + C2 * log(I₂(λ⃗) / 3)
end
# Suitable for low strains
"""
Hoss Marczak I

Parameters: α, β, μ, b, n

Model: ``\\frac{\\alpha}{\\beta}(1-\\exp{-\\beta(I_1-3)})+\\frac{\\mu}{2b}\\bigg((1+\\frac{b}{n}(I_1-3))^n -1\\bigg)``
"""
function HossMarczakI((; α, β, μ, b, n))
    (λ⃗) -> α / β * (1 - exp(-β * (I₁(λ⃗) - 3))) + μ / (2b) * ((1 + b / n * (I₁(λ⃗) - 3))^n - 1)
end
# Suitable for high strains
"""
Hoss Marczak II

Parameters: α, β, μ, b, n, C2

Model: ``\\frac{\\alpha}{\\beta}(1-\\exp{-\\beta(I_1-3)})+\\frac{\\mu}{2b}\\bigg((1+\\frac{b}{n}(I_1-3))^n -1\\bigg)+C_2\\log(\\frac{I_2}{3})``
"""
function HossMarczakII((; α, β, μ, b, n, C2))
    (λ⃗) -> α / β * (1 - exp(-β * (I₁(λ⃗) - 3))) + μ / (2b) * ((1 + b / n * (I₁(λ⃗) - 3))^n - 1) + C2 * log(I₂(λ⃗) / 3)
end

"""
Exp-Ln

Parameters: A, a, b

Model: ``A\\bigg[\\frac{1}{a}\\exp{(a(I_1-3))}+b(I_1-2)(1-\\log{I_1-2})-\\frac{1}{a}-b\\bigg]``
"""
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

"""
Gent

Parameters: μ, Jₘ

Model: ``-\\frac{\\mu J_m}{2}\\log{1-\\frac{I_1-3}{J_m}}``
"""
function Gent((; μ, Jₘ))
    (λ⃗) -> -μ * Jₘ / 2 * log(1 - (I₁(λ⃗) - 3) / Jₘ)
end

# With the assumption of isotropicity -> Verified with A description of arterial wall mechanics using limiting chain extensibility constitutitive models by Horgan and Saccomandi
"""
Takamizawa-Hayashi

Parameters: c, Jₘ

Model: ``-c\\log{1-\\big(\\frac{I_1-3}{J_m}\\big)^2}``
"""
function TakamizawaHayashi((; c, Jₘ))
    (λ⃗) -> -c * log(1 - ((I₁(λ⃗) - 3) / Jₘ)^2)
end

function YeohFleming((; A, B, C10, Im))
    (λ⃗) -> A / B * (1 - exp(-B * (I₁(λ⃗) - 3))) - C10 * (Im - 3) * log(1 - ((I₁(λ⃗) - 3) / (Im - 3)))
end

## Not a real model => Not referenced in the Gent paper cited
# function Gent3Parameters((;μ, Jₘ, α))
#     (λ⃗) -> 
# end

"""
Pucci-Saccomandi

Parameters: K, μ, Jₘ

Model ``K\\log{\\frac{I_2}{3}}-\\frac{\\mu J_m}{2}\\log{1-\\frac{I_1-3}{J-m}}``
"""
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
########### TABLE 4
########################
"""
Valanis-Landel

Parameters: μ

Model: ``2\\mu\\sum\\limits_{1}^{3}(\\lambda_i(\\log\\lambda_i -1))``
"""
function ValanisLandel((; μ))
    (λ⃗) -> 2 * μ * sum(λ⃗ * (log.(λ⃗) - 1))
end

"""
Peng - Landel

Parameters: E

Model: ``E\\sum\\limits_{i=1}^{3}\\bigg[\\lambda_i - 1 - \\log(\\lambda_i) - \\frac{1}{6}\\log(\\lambda_i)^2 + \\frac{1}{18}\\log(\\lambda_i)^3-\\frac{1}{216}\\log(\\lambda_i)^4]``
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

Model: ``\\sum\\limits_{i=1}^{N} A_i\\big[\\exp{m_i(\\lambda_1^{\\alpha_i}+\\lambda_2^{\\alpha_i}+\\lambda_3^{\\alpha_i}-3)}-1] + B_i\\big[\\exp{n_i(\\lambda_1^{-\\beta_i}+\\lambda_2^{-\\beta_i}+\\lambda_3^{-\\beta_i}-3)}-1]``
"""
function ArmanNarooei((; A, B, m, n, α, β))
    (λ⃗) -> @tullio _ := A[i] * (exp(m[i] * (sum(λ⃗ .^ α[i]) - 3)) - 1) + B[i] * (exp(n[i] * (sum(λ⃗ .^ (-β[i])) - 3)) - 1)
end

################################
###################      Table 5
################################
"""
Continuum Hybrid

Parameters: K₁, K₂, α, μ

Model: ``K_1(I_1-3)+K_2\\log\\frac{I_2}{3}+\\frac{\\mu}{\\alpha}(\\lambda_1^\\alpha+\\lambda_2^\\alpha+\\lambda^\\alpha-3)``
"""
function ContinuumHybrid((; K₁, K₂, α, μ))
    (λ⃗) -> K₁ * (I₁(λ⃗) - 3) + K₂ * log(I₂ / 3) + μ / α * (sum(λ⃗ .^ α) - 3)
end

"""
Bechir-4 Term

Parameters: C11, C12, C21, C22

Model: ``C_1^1(I_1-3)+\\sum\\limits_{n=1}^{2}\\sum\\limits_{r=1}^{2}C_n^{r}(\\lambda_1^{2n}+\\lambda_2^{2n}+\\lambda_3^{2n}-3)^r``
"""
function Bechir4Term((; C11, C12, C21, C22))
    C = [C11 C12; C21 C22]
    (λ⃗) -> C[1, 1] * (I₁(λ⃗) - 3) + sum(n -> sum(r -> C[n, r] * (sum(λ⃗ .^ (2n))), 1:2), 1:2)
end

"""
WFB - Skipped

Parameters: Lf, F, A, B, C, D

Model: ``\\int\\limits_{1}^{L_f}\\big(F(\\lambda_1)A()\\big)``
"""
function WFB()
    error("Not Yet Implemented")
end

"""
Constrained Junction

Parameters: Gc, νkT, κ  

Model:
"""
function ConstrainedJunction((; Gc, νkT, κ))
    (λ⃗) -> Gc * (I₁(λ⃗) - 3) + μkT / 2 * sum(i -> κ * (λ⃗[i] - 1) / (λ⃗[i]^2 + κ) + log((λ⃗[i]^2 + κ) / (1 + κ)) - log(λ⃗[i]^2), 1:3)
end
"""
Edward-Vilgis

Parameters: Ns, Nc, α, η

Model: ``\\frac{1}{2}N_C\\Bigg[\\frac{(1-\\alpha^2)I_1}{1-\\alpha^2I_1}+\\log(1-\\alpha^2I_1)\\Bigg]+\\frac{1}{2}N_S\\Bigg[\\sum_{i=1}^{3}\\Big\\{\\frac{(1+\\eta)(1-\\alpha^2)\\lambda_i^2}{( 1+\\eta\\lambda_i^2)(1-\\alpha^2I_1)}+\\log(1+\\eta\\lambda_i^2)\\Big\\}+\\log(1-\\alpha^2I_1)\\Bigg]``
"""
function EdwardVilgis((; Ns, Nc, α, η))
    function W(λ⃗)
        A = 0.5 * Nc * ((1 - α^2) * I₁(λ⃗) / (1 - α^2 * I₁(λ⃗)) + log(1 - α^2 * I₁(λ⃗)))
        B = sum(i -> (1 + η) * (1 - α^2) * λ⃗[i] / (1 - η * λ⃗[i]^2) / (1 - α^2 * I₁(λ⃗)) + log(1 + η * λ⃗[i]^2), 1:3)
        B = 0.5 * Ns * (B + log(1 - α^2 * I₁(λ⃗)))
        W = A + B
        return W
    end
end

"""
MCC (modified constrained chain)

Parameters:

Model:``\\frac{1}{2}\\zeta k T \\sum\\limits_{i=1}^{3}(\\lambda_i^2-1)+\\frac{1}{2}\\mu k T\\sum\\limits_{i=1}^{3}[B_i+D_i-\\log{(1+B_i)}-\\log{(1+D_i)}]``
``B_i = \\frac{\\kappa^2(\\lambda_i^2-1)}{(\\lambda_i^2+\\kappa)^2}``
``D_i = \\frac{\\lambda_i^2 B_i}{\\kappa}``
"""
function MCC((; ζkT, μkT, κ))
    (λ⃗) ->
        1 / 2 * ζkT * sum(i -> λ⃗[i]^2 - 1, 1:3) + 1 / 2 * μkT * sum(i -> κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2) + (λ⃗[i]^2 * (κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2)) / κ) - log(1 + (κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2))) - log(1 + (λ⃗[i]^2 * (κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2)) / κ)))
end

"""
Tube

Parameters: Gc, Ge, β

Model: ``\\sum\\limits_{i=1}^{3}\\frac{G_c}{2}(\\lambda_i^2-1)+\\frac{2Ge}{\\beta^2}(\\lambda_i^{-\\beta}-1)``
"""
function Tube((; Gc, Ge, β))
    (λ⃗) -> @tullio _ := Gc / 2 * (λ⃗[i]^2 - 1) + 2Ge / β^2 * (λ⃗[i]^(-β) - 1)
end

"""
Nonaffine - Tube

Parameters: Gc, Ge

Model: ``G_c \\sum\\limits_{i=1}^{3}\\frac{\\lambda_i^2}{2}+G_e\\sum\\limits_{i=1}^{3}\\lambda_i+\\frac{1}{\\lambda_i}``
"""
function NonaffineTube((; Gc, Ge))
    (λ⃗) -> Gc * sum(λ⃗ .^ 2 ./ 2) + Ge * sum(λ⃗ .+ 1 ./ λ⃗)
end
# ---------------------------------------------------- #

"""
Arruda Boyce 

Parameters: μ, N

Model: ``\\mu\\bigg(\\frac{1}{2}(I_1-3)+\\frac{I_1^2-9}{20N}+\\frac{11(I_1^3-27)}{1050N^2}+\\frac{19(I_1^4-81)}{7000N^3}+\\frac{519(I_1^5-243)}{673750N^4}\\bigg)``
"""
function ArrudaBoyce((; μ, N))
    (λ⃗) -> μ * (0.5 * (I₁(λ⃗) - 3) + 1 / 20 / N * (I₁(λ⃗)^2 - 9) + 11 / 1050 / N^2 * (I₁(λ⃗) - 27) + 19 / 7000 / N^3 * (I₁(λ⃗)^4 - 81) + 519 / 673750 / N^4 * (I₁(λ⃗)^5 - 243))
end

end