# # Available Models
export GeneralMooneyRivlin, GeneralDarijaniNaghdabadi, GeneralBeda, MooneyRivlin, NeoHookean, Gent, Biderman, Isihara, JamesGreenSimpson, Lion, Yeoh, HauptSedlan, HartmannNeff, HainesWilson, Carroll, BahremanDarijani, Zhao, Knowles, Swanson, YamashitaKawabata, DavisDeThomas, Gregory, ModifiedGregory, Beda, Amin, LopezPamies, GenYeoh, HartSmith, VerondaWestmann, FungDemiray, Vito, ModifiedYeoh, Martins, ChevalierMarco, GornetDesmorat, MansouriDarijani, GentThomas, Alexander, LambertDianiRey, HossMarczakI, HossMarczakII, ExpLn, Kilian, VanDerWaals, TakamizawaHayashi, YeohFleming, PucciSaccomandi, HorganSaccomandi, Beatty, HorganMurphy, ArrudaBoyce, Ogden, EdwardVilgis, NonaffineTube, Tube, MCC, Bechir4Term, ConstrainedJunction, ContinuumHybrid, ArmanNarooei, PengLandel, ValanisLandel, Attard, Shariff, ThreeChainModel, ModifiedFloryErman, ABGI, BechirChevalier, Bootstrapped8Chain, DavidsonGoulbourne, ExtendedTubeModel, FullNetwork, GeneralConstitutiveModel, Lim, MicroSphere, NetworkAveragingTube, WFB, ZunigaBeatty

"""
General Mooney Rivlin

Parameters: [C]

Model:    
``\\sum\\limits_{i,j = 0}^{N,M} C_{i,j}(I_1-3)^i(I_2-3)^j``
"""
function GeneralMooneyRivlin((; C))
    function W(λ⃗)
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
    W(λ⃗) = sum(A⃗ .* (λ⃗ .^ m⃗ .- 1) + B⃗ .* (λ⃗ .^ (-n⃗) .- 1))
end

"""
General Beda

Parameters: C, K, α, β

Model: ``\\sum\\limits_{i = 1}^{N}\\frac{C_i}{\\alpha_i}(I_1-3)^{\\alpha_i} + \\sum\\limits_{j=1}^{M}\\frac{K_j}{\\beta_j}(I_2-3)^{\\beta_j}``
"""
function GeneralBeda((; C, K, α, β))
    function W(λ⃗)
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
end

"""
NeoHookean

Parameters: μ

Model: ``\\frac{\\mu}{2}(I_1-3)``
"""
function NeoHookean((; μ))
    W = GeneralMooneyRivlin(ComponentVector(
        C=[
        0.0 μ
    ]))
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
end

"""
Haupt Sedlan

Parameters: C10, C01, C11, C02, C30

Model: 
``\\sum\\limits_{i,j=0}^{3, 2}C_{i,j}(I_1-3)^i(I_2-3)^j``
"""
function HauptSedlan((; C10, C01, C11, C02, C30))
    W = GeneralMooneyRivlin(ComponentVector(
        C=[
            0.0 C10 0.0 C30
            C01 C11 0.0 0.0
            C02 0.0 0.0 0.0
        ]))
end

"""
Hartmann-Neff

Parameters: α, Ci0, C0j

Model: ``\\sum\\limits_{i,j=0}^{M,N}C_{i,0}(I_1-3)^i -3\\sqrt{3}^j+\\alpha(I_1-3)``
"""
function HartmannNeff((; α, Ci0, C0j))
    function W(λ⃗)
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
    W(λ⃗) = A * I₁(λ⃗) + B * I₁(λ⃗)^4 + C * I₂(λ⃗)^(1 / 2)
end

## Only developed for simple shear deformation.
# function Nunes((; C1, C2))
# W(λ⃗) = C1 * (I₁(λ⃗) - 3) + 4 / 3 * C2 * (I₂(λ⃗) - 3)^(3 / 4)
# W(λ⃗) = 1/2*(C1*(I₁(λ⃗)-3)+C2*(I₂(λ⃗)-3)^(3/4))
# end

function BahremanDarijani((; A2, B2, A4, A6))
    W = GeneralizedDarijaniNaghdabadi(
        ComponentVector(
            A=[0.0, A2, 0.0, A4, 0.0, A6],
            B=[0.0, B2],
            m=[0.0, 2.0, 0.0, 4.0, 0.0, 6.0],
            n=[0.0, 2.0])
    )
end
"""
Zhao

Parameters: C₋₁¹,, C₁¹, C₂¹, C₂²

Model: ``C_{-1}^1*(I_2-3)+C_{1}^{1}(I_1-3)+C_{2}^{1}(I_1^2-2I_2-3)+C_{2}^{2}(I_1^2-2I_2-3)^2``
"""
function Zhao((; C₋₁¹, C₁¹, C₂¹, C₂²))
    W(λ⃗) = C₋₁¹ * (I₂(λ⃗) - 3) + C₁¹ * (I₁(λ⃗) - 3) + C₂¹ * (I₁(λ⃗)^2 - 2I₂(λ⃗) - 3) + C₂² * (I₁(λ⃗)^2 - 2I₂(λ⃗) - 3)^2
end

## Table 2
"""
Knowles

Parameters: μ, b, n

Model: ``\\frac{\\mu}{2b}((1+\\frac{b}{n}(I_1-3))^n-1)``
"""
function Knowles((; μ, b, n))
    W(λ⃗) = μ / (2b) * ((1 + (b / n) * (I₁(λ⃗) - 3))^n - 1)
end

# Article Requested
"""
Swanson

Parameters: A, α, B, β

Model: ``\\sum\\limits_{i=1}^{N} \\frac{3}{2}(\\frac{A_i}{1+\\alpha_i}(\\frac{I_1}{3})^{1+\\alpha_i}+\\frac{B_i}{1+\\beta_i}(\\frac{I_2}{3})^{1+\\beta_i}``
"""
function Swanson((; A, α, B, β))
    W(λ⃗) = @tullio _ := 3 / 2 * (A[i] / (1 + α[i]) * (I₁(λ⃗) / 3)^(1 + α[i]) + B[i] / (1 + β[i]) * (I₂(λ⃗) / 3)^(1 + β[i]))
end

# Original article in Japanese
"""
Yamashita-Kawabata

Parameters: C1, C2, C3, N

Model: ``C_1(I_1-3)+C_2(I_2-3)+\\frac{C_3}{N+1}(I_1-3)^{N+1}``
"""
function YamashitaKawabata((; C1, C2, C3, N))
    W(λ⃗) = C1 * (I₁(λ⃗) - 3) + C2 * (I₂(λ⃗) - 3) + C3 / (N + 1) * (I₁(λ⃗) - 3)^(N + 1)
end

# Article Requested
"""
Davis-DeThomas

Parameters: A, n, C, k

Model: ``\\frac{A}{2(1-\\frac{n}{2})}(I_1-3+C^2)^{1-\\frac{n}{2}}+k(I_1-3)^2``
"""
function DavisDeThomas((; A, n, C, k))
    W(λ⃗) = A / (2 * (1 - n / 2)) * (I₁(λ⃗) - 3 + C^2)^(1 - n / 2) + k * (I₁(λ⃗) - 3)^2
end

# Article Requested
"""
Gregory

Parameters: A, B, C, m, n

Model: ``\\frac{A}{2-n}(I_1-3+C^2)^{1-\\frac{n}{2}}+\\frac{B}{2+m}(I_1-3+C^2)^{1+\\frac{m}{2}}``
"""
function Gregory((; A, B, C, m, n))
    W(λ⃗) = A / (2 - n) * (I₁(λ⃗) - 3 + C^2)^(1 - n / 2) + B / (2 + m) * (I₁(λ⃗) - 3 + C^2)^(1 + m / 2)
end

# Proposed in 85 Model review
"""
Modified Gregory

Parameters: A, α, M, B, β, N

Model: ``\\frac{A}{1+\\alpha}(I_1-3+M^2)^{1+\\alpha}+\\frac{B}{1+\\beta}(I_1-3+N^2)^{1+\\beta}``
"""
function ModifiedGregory((; A, α, M, B, β, N))
    W(λ⃗) = A / (1 + α) * (I₁(λ⃗) - 3 + M^2)^(1 + α) + B / (1 + β) * (I₁(λ⃗) - 3 + N^2)^(1 + β)
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
    # W(λ⃗) = C1 / α * (I₁(λ⃗) - 3)^(α) + C2 * (I₁(λ⃗) - 3) + C3 / ζ * (I₁(λ⃗) - 3)^(ζ) + K1 / β * (I₂(λ⃗) - 3)^β
end

"""
Amin

Parameters: C1, C2, C3, C4, N, M

Model:``C_1 (I_1 - 3) + \\frac{C_2}{N + 1} (I_1 - 3)^{N + 1} + \\frac{C_3}{M + 1} (I_1 - 3)^{M + 1} + C_4 (I_2 - 3)``
"""
function Amin((; C1, C2, C3, C4, N, M))
    W(λ⃗) = C1 * (I₁(λ⃗) - 3) + C2 / (N + 1) * (I₁(λ⃗) - 3)^(N + 1) + C3 / (M + 1) * (I₁(λ⃗) - 3)^(M + 1) + C4 * (I₂(λ⃗) - 3)
end

# Article Received - General form presented.
"""
Lopez-Pamies

Parameters: α⃗, μ⃗

Model: ``\\frac{3.0^{1 - \\alpha_i}}{2\\alpha_i} \\mu_i (I_1^{\\alpha_i} - 3^{\\alpha_i})``
"""
function LopezPamies((; α⃗, μ⃗))
    W(λ⃗) = @tullio _ := (3.0^(1 - α⃗[i])) / (2α⃗[i]) * μ⃗[i] * (I₁(λ⃗)^(α⃗[i]) - 3^(α⃗[i]))
end

# ✔
"""
GenYeoh

Parameters: K1, K2, K3, m, p, q

Model: ``K_1 (I_1 - 3)^m + K_2 * (I_1 - 3)^p + K_3 * (I_1 - 3)^q``
"""
function GenYeoh((; K1, K2, K3, m, p, q))
    W(λ⃗) = K1 * (I₁(λ⃗) - 3)^m + K2 * (I₁(λ⃗) - 3)^p + K3 * (I₁(λ⃗) - 3)^q
end
# ✔
"""
Veronda-Westmann

Parameters: C1, C2, α

Model: ``C_1 (\\exp(\\alpha(I_1 - 3)) - 1) + C_2 (I_2 - 3)``
"""
function VerondaWestmann((; C1, C2, α))
    W(λ⃗) = C1 * (exp(α * (I₁(λ⃗) - 3)) - 1) + C2 * (I₂(λ⃗) - 3)
end
# ✔
"""
Fung-Demiray

Parameters: μ, b

Model: ``\\frac{\\mu}{2 * b} (\\exp(b(I_1 - 3)) - 1)``
"""
function FungDemiray((; μ, b))
    W(λ⃗) = μ / (2 * b) * (exp(b * (I₁(λ⃗) - 3)) - 1)
end
# ✔
"""
Vito

Parameters: α, β, γ

Model: ``\\alpha (\\exp\\bigg(\\beta (I_1 - 3)\\bigg) + \\gamma  (I_2 - 3)) - 1)``
"""
function Vito((; α, β, γ))
    W(λ⃗) = α * (exp(β * (I₁(λ⃗) - 3) + γ * (I₂(λ⃗) - 3)) - 1)
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
    W(λ⃗) = C10 * (I₁(λ⃗) - 3) + C20 * (I₁(λ⃗) - 3)^2 + C30 * (I₁(λ⃗) - 3)^3 + α / β * (1 - exp(-β * (I₁(λ⃗) - 3)))
end
# ✔
"""
Mansouri-Darijani

Parameters: A1, m1, B1, n1

Model: ``A_1\\exp{m_1(I_1-3)-1}+B_1\\exp{n_1(I_2-3)-1}``
"""
function MansouriDarijani((; A1, m1, B1, n1))
    W(λ⃗) = A1 * (exp(m1 * (I₁(λ⃗) - 3)) - 1) + B1 * (exp(n1 * (I₂(λ⃗) - 3)) - 1)
end
# ✔
"""
Gent Thomas

Paramters: C1, C2

Model: ``C_1(I_1-3)+C_2\\log(\\frac{I_2}{3})``
"""
function GentThomas((; C1, C2))
    W(λ⃗) = C1 * (I₁(λ⃗) - 3) + C2 * log(I₂(λ⃗) / 3)
end
# Suitable for low strains
"""
Hoss Marczak I

Parameters: α, β, μ, b, n

Model: ``\\frac{\\alpha}{\\beta}(1-\\exp{-\\beta(I_1-3)})+\\frac{\\mu}{2b}\\bigg((1+\\frac{b}{n}(I_1-3))^n -1\\bigg)``
"""
function HossMarczakI((; α, β, μ, b, n))
    W(λ⃗) = α / β * (1 - exp(-β * (I₁(λ⃗) - 3))) + μ / (2b) * ((1 + b / n * (I₁(λ⃗) - 3))^n - 1)
end
# Suitable for high strains
"""
Hoss Marczak II

Parameters: α, β, μ, b, n, C2

Model: ``\\frac{\\alpha}{\\beta}(1-\\exp{-\\beta(I_1-3)})+\\frac{\\mu}{2b}\\bigg((1+\\frac{b}{n}(I_1-3))^n -1\\bigg)+C_2\\log(\\frac{I_2}{3})``
"""
function HossMarczakII((; α, β, μ, b, n, C2))
    W(λ⃗) = α / β * (1 - exp(-β * (I₁(λ⃗) - 3))) + μ / (2b) * ((1 + b / n * (I₁(λ⃗) - 3))^n - 1) + C2 * log(I₂(λ⃗) / 3)
end

"""
Exp-Ln

Parameters: A, a, b

Model: ``A\\bigg[\\frac{1}{a}\\exp{(a(I_1-3))}+b(I_1-2)(1-\\log{I_1-2})-\\frac{1}{a}-b\\bigg]``
"""
function ExpLn((; A, a, b))
    W(λ⃗) = A * (1 / a * exp(a * (I₁(λ⃗) - 3)) + b * (I₁(λ⃗) - 2) * (1 - log(I₁(λ⃗) - 2)) - 1 / a - b)
end

#####################
###### TABLE 3 ######
#####################
# model not mentioned in original article ?????? 
# function Warner((; μ, Iₘ))
#     W(λ⃗) = -1/2*μ*Iₘlog(1-(I₁(λ⃗)-3)/(Iₘ-3))
# end 

# does not easily match the form in the paper
# function Killian((; μ, JL))
#     W(λ⃗) = -μ * JL * (log(1 - sqrt((I₁(λ⃗) - 3) / JL)) + sqrt((I₁(λ⃗) - 3) / JL))
# end
# ARticles requested -> Checked against review article from Marckmann and Verron
function VanDerWaals((; μ, λm, β, α))
    W(λ⃗) = μ * (-(λm^2 - 3) * (log(1 - θ) + θ) - 2 / 3 * α * ((I₁(λ⃗) - 3) / 2)^(3 / 2))
end

"""
Gent

Parameters: μ, Jₘ

Model: ``-\\frac{\\mu J_m}{2}\\log{\\bigg(1-\\frac{I_1-3}{J_m}\\bigg)}``
"""
function Gent((; μ, Jₘ))
    W(λ⃗) = -(μ * Jₘ) / 2 * log(1 - (I₁(λ⃗) - 3) / Jₘ)
end

"""
Takamizawa-Hayashi
From: A description of arterial wall mechanics using limiting chain extensibility constitutitive models by Horgan and Saccomandi
Parameters: c, Jₘ

Model: ``-c\\log{1-\\big(\\frac{I_1-3}{J_m}\\big)^2}``
"""
function TakamizawaHayashi((; c, Jₘ))
    W(λ⃗) = -c * log(1 - ((I₁(λ⃗) - 3) / Jₘ)^2)
end

function YeohFleming((; A, B, C10, Im))
    W(λ⃗) = A / B * (1 - exp(-B * (I₁(λ⃗) - 3))) - C10 * (Im - 3) * log(1 - ((I₁(λ⃗) - 3) / (Im - 3)))
end

## Not a real model => Not referenced in the Gent paper cited
# function Gent3Parameters((;μ, Jₘ, α))
#     W(λ⃗) = 
# end

"""
Pucci-Saccomandi

Parameters: K, μ, Jₘ

Model ``K\\log{\\frac{I_2}{3}}-\\frac{\\mu J_m}{2}\\log{1-\\frac{I_1-3}{J-m}}``
"""
function PucciSaccomandi((; K, μ, Jₘ))
    W(λ⃗) = K * log(I₂(λ⃗) / 3) - μ * Jₘ / 2 * log(1 - (I₁(λ⃗) - 3) / Jₘ)
end

# Originally from CONSTITUTIVE MODELS FOR ATACTIC ELASTOMERS
"""
Horgan Saccomandi Model

Parameters: μ, J

Model: ``-\\frac{\\mu J}{2}\\log\\bigg(\\frac{J^3-J^2I_1+JI_2-1}{(J-1)^3}\\bigg)``
"""
function HorganSaccomandi((; μ, J))
    W(λ⃗) = -μ * J / 2 * log((J^3 - J^2 * I₁(λ⃗) + J * I₂(λ⃗) - 1) / (J - 1)^3)
end

"""
Beatty Model

Parameters: G₀, Iₘ

Model: ``-\\frac{G_0 I_m(I_m-3)}{2(2I_m-3)}\\log\\bigg(\\frac{1-\\frac{I_1-3}{I_m-3}}{1+\\frac{I_1-3}{I_m}} \\bigg)``
"""
function Beatty((; G₀, Iₘ))
    W(λ⃗) = -G₀ * Iₘ * (Iₘ - 3) / 2 / (2Iₘ - 3) * log((1 - (I₁(λ⃗) - 3) / (Iₘ - 3)) / (1 + (I₁(λ⃗) - 3) / (Iₘ)))
end

"""
Horgan Murphy Model

Parameters: μ, Jₘ, c

Model: ``-\\frac{2\\mu J_m}{c^2}\\log\\bigg(1-\\frac{\\lambda_1^c+\\lambda_2^c+\\lambda_3^c-3}{J_m})``
"""
function HorganMurphy((; μ, Jₘ, c))
    W(λ⃗) = -2 * μ * Jₘ / c^2 * log(1 - (sum(λ⃗ .^ c) - 3) / Jₘ)
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
    W(λ⃗) = 2 * μ * sum(λ⃗ * (log.(λ⃗) - 1))
end

"""
Peng - Landel

Parameters: E

Model: ``E\\sum\\limits_{i=1}^{3}\\bigg[\\lambda_i - 1 - \\log(\\lambda_i) - \\frac{1}{6}\\log(\\lambda_i)^2 + \\frac{1}{18}\\log(\\lambda_i)^3-\\frac{1}{216}\\log(\\lambda_i)^4]``
"""
function PengLandel((; E))
    W(λ⃗) = sum(@. λ⃗ - 1 - log(λ⃗) - 1 / 6 * log(λ⃗)^2 + 1 / 18 * log(λ⃗)^3 - 1 / 216 * log(λ⃗)^4) * E
end

"""
Ogden

Parameters: μ⃗, α⃗

Model: ``\\sum\\limits_{i=1}^{N}\\frac{\\mu_i}{\\alpha_i}(\\lambda_1^{\\alpha_i}+\\lambda_2^{\\alpha_i}+\\lambda_3^{\\alpha_i}-3)``
"""
function Ogden((; μ, α))
    W(λ⃗) = @tullio _ := μ[i] / α[i] * (sum(λ⃗ .^ α[i]) - 3)
end

"""
Attard

Parameters: A⃗, B⃗

Model: ``\\sum\\limits_{i=1}^N\\frac{A_i}{2i}(\\lambda_1^{2i}+\\lambda_2^{2i}+\\lambda_3^{2i}-3) + \\frac{B_i}{2i}(\\lambda_1^{-2i}+\\lambda_2^{-2i}+\\lambda_3^{-2i}-3)``
"""
function Attard((; A, B))
    W(λ⃗) = @tullio _ := A[i] / 2 / i * (sum(λ⃗ .^ (2i)) - 3) + B[i] / 2 / i * (sum(λ⃗ .^ (-2i)) - 3)
end

"""
Shariff

Parameters: E, α₁, α₂, α₃, α₄, α₅

Model: 
``E\\sum\\limits_{i=1}^3\\sum\\limits_{j=1}^{N}\\alpha_j \\Phi_j(\\lambda_i)``
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
    W(λ⃗) = E * (@tullio _ := ϕ[i](λ⃗[j]))
end

# Article requested
"""
Arman - Narooei

Parameters: A⃗, B⃗, m⃗, n⃗, α⃗, β⃗

Model: ``\\sum\\limits_{i=1}^{N} A_i\\big[\\exp{m_i(\\lambda_1^{\\alpha_i}+\\lambda_2^{\\alpha_i}+\\lambda_3^{\\alpha_i}-3)}-1] + B_i\\big[\\exp{n_i(\\lambda_1^{-\\beta_i}+\\lambda_2^{-\\beta_i}+\\lambda_3^{-\\beta_i}-3)}-1]``
"""
function ArmanNarooei((; A, B, m, n, α, β))
    W(λ⃗) = @tullio _ := A[i] * (exp(m[i] * (sum(λ⃗ .^ α[i]) - 3)) - 1) + B[i] * (exp(n[i] * (sum(λ⃗ .^ (-β[i])) - 3)) - 1)
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
    W(λ⃗) = K₁ * (I₁(λ⃗) - 3) + K₂ * log(I₂ / 3) + μ / α * (sum(λ⃗ .^ α) - 3)
end

"""
Bechir-4 Term

Parameters: C11, C12, C21, C22

Model: ``C_1^1(I_1-3)+\\sum\\limits_{n=1}^{2}\\sum\\limits_{r=1}^{2}C_n^{r}(\\lambda_1^{2n}+\\lambda_2^{2n}+\\lambda_3^{2n}-3)^r``
"""
function Bechir4Term((; C11, C12, C21, C22))
    C = [C11 C12; C21 C22]
    W(λ⃗) = C[1, 1] * (I₁(λ⃗) - 3) + sum(n -> sum(r -> C[n, r] * (sum(λ⃗ .^ (2n))), 1:2), 1:2)
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

Model: ``G_c (I_1-3)+ \\frac{\\nu k T}{2}(\\sum\\limits_{i=1}^{3}\\kappa\\frac{\\lambda_i-1}{\\lambda_i^2+\\kappa}+\\log{\\frac{\\lambda_i^2+\\kappa}{1+\\kappa}}-\\log{\\lambda_i^2})``
"""
function ConstrainedJunction((; Gc, νkT, κ))
    W(λ⃗) = Gc * (I₁(λ⃗) - 3) + μkT / 2 * sum(i -> κ * (λ⃗[i] - 1) / (λ⃗[i]^2 + κ) + log((λ⃗[i]^2 + κ) / (1 + κ)) - log(λ⃗[i]^2), 1:3)
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
    W(λ⃗) =
        1 / 2 * ζkT * sum(i -> λ⃗[i]^2 - 1, 1:3) + 1 / 2 * μkT * sum(i -> κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2) + (λ⃗[i]^2 * (κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2)) / κ) - log(1 + (κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2))) - log(1 + (λ⃗[i]^2 * (κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2)) / κ)))
end

"""
Tube

Parameters: Gc, Ge, β

Model: ``\\sum\\limits_{i=1}^{3}\\frac{G_c}{2}(\\lambda_i^2-1)+\\frac{2Ge}{\\beta^2}(\\lambda_i^{-\\beta}-1)``
"""
function Tube((; Gc, Ge, β))
    W(λ⃗) = @tullio _ := Gc / 2 * (λ⃗[i]^2 - 1) + 2Ge / β^2 * (λ⃗[i]^(-β) - 1)
end

"""
Nonaffine - Tube

Parameters: Gc, Ge

Model: ``G_c \\sum\\limits_{i=1}^{3}\\frac{\\lambda_i^2}{2}+G_e\\sum\\limits_{i=1}^{3}\\lambda_i+\\frac{1}{\\lambda_i}``
"""
function NonaffineTube((; Gc, Ge))
    W(λ⃗) = Gc * sum(λ⃗ .^ 2 ./ 2) + Ge * sum(λ⃗ .+ 1 ./ λ⃗)
end
"""
Three Chain Model

Parameters: μ, N

Model: `` \\frac{\\mu\\sqrt{N}}{3}\\sum\\limits_{i=1}^{3}\\bigg(\\lambda_i\\beta_i+\\sqrt{N}\\log\\bigg(\\frac{\\beta_i}{\\sinh \\beta_i}\\bigg)\\bigg)

"""
function ThreeChainModel((; μ, N))
    ℒinv(x) = x * (3 - 1.0651 * x^2 - 0.962245 * x^4 + 1.47353 * x^6 - 0.48953 * x^8) / (1 - x) / (1 + 1.01524 * x)
    W(λ⃗) = μ * sqrt(N) / 3 * sum(λ⃗ .* ℒinv.(λ⃗ ./ sqrt(N)) .+ sqrt(N) .* log.((ℒinv.(λ⃗ ./ sqrt(N))) ./ (sinh.(ℒinv.(λ⃗ ./ sqrt(N))))))
end

"""
Arruda Boyce 

Parameters: μ, N

Model: ``\\mu\\bigg(\\frac{1}{2}(I_1-3)+\\frac{I_1^2-9}{20N}+\\frac{11(I_1^3-27)}{1050N^2}+\\frac{19(I_1^4-81)}{7000N^3}+\\frac{519(I_1^5-243)}{673750N^4}\\bigg)``
"""
function ArrudaBoyce((; μ, N))
    ℒinv(x) = x * (3 - 1.0651 * x^2 - 0.962245 * x^4 + 1.47353 * x^6 - 0.48953 * x^8) / (1 - x) / (1 + 1.01524 * x)
    function W(λ⃗)
        rchain_Nl = √(I₁(λ⃗) / 3 / N)
        β = ℒinv(rchain_Nl)
        μ * N * (rchain_Nl * β + log(β / sinh(β)))
    end
end

"""
Modified Flory Erman

Parameters: μ, N, κ

Model: ``W_{\\text{Arruda-Boyce}}+\\sum\\limits_{i=1}^{3}\\frac{\\mu}{2}[B_i+D_i]
"""

function ModifiedFloryErman((; μ, N, κ))
    function W(λ⃗)
        B = map(i -> κ^2 * (λ⃗[i]^2 - 1) / (λ⃗[i]^2 + κ)^2, 1:3)
        D = map(i -> λ⃗[i]^2 * B[i] / κ, 1:3)
        ArrudaBoyce((μ=μ, N=N))(λ⃗) + map(i -> B[i] + D[i] - log(B[i] + 1) - log(D[i] + 1), 1:3)
    end
end

"""
Extended Tube Model

Parameters: Gc, Ge, δ, β

Model: ``\\frac{G_c}{2}\\bigg[\\frac{(1-\\delta^2)(I_1-3)}{1-\\delta^2(I_1-3)}+\\log{(1-\\delta^2(I_1-3))}\\bigg]+\\frac{2G_e}{\\beta^2}\\sum\\limits_{i=1}^{3}(\\lambda_i^{-\\beta}-1)
"""
function ExtendedTubeModel((Gc, Ge, δ, β))
    W(λ⃗) = Gc / 2 * ((1 - δ^2) * (I₁(λ⃗) - 3) / (1 - δ^2 * (I₁ - 3)) + log(1 - δ^2 * (I₁(λ⃗) - 3))) + 2 * Ge / β^2 * sum(λ⃗ .^ (-β) .- 1)
end

"""
ABGI

Parameters: μ, N, Ge, n

Model: ``W_{Arruda-Boyce} + G_e\\frac{\\lambda_1^n+\\lambda_2^2+\\lambda_3^2-3}{n}``
"""
function ABGI((; μ, N, Ge, n))
    W(λ⃗) = ArrudaBoyce((μ=μ, N=N))(λ⃗) + Ge * (sum(λ⃗ .^ n) - 3) / n
end

"""
Micro-Sphere

Parameters:

Model:
"""
function MicroSphere((; μ, N, P, U, q))
    error("Not Yet implemented")
end


"""
Bootstrapped 8Chain Model

Parameters: μ, N

Model: ``W_8(\\frac{\\sum\\lambda}{\\sqrt{3N}}-\\frac{\\lambda_{chain}}{\\sqrt{N}})+W_{8}(\\frac{\\lambda_{chain}}{\\sqrt{N}})``

``W_8(x) = \\mu N (x \\mathcal{L}^{-1}(x) + \\log\\frac{\\mathcal{L}^{-1}(x)}{\\sinh\\mathcal{L}^{-1}(x)})``

``\\lambda_{cha`in} = \\sqrt{\\frac{I_1}{3}}``
"""
function Bootstrapped8Chain((; μ, N))
    ℒinv(x) = x * (3 - 1.0651 * x^2 - 0.962245 * x^4 + 1.47353 * x^6 - 0.48953 * x^8) / (1 - x) / (1 + 1.01524 * x)
    W8(x) = μ * N * (x * ℒinv(x) + log(ℒinv(x) / sinh(ℒinv(x))))
    function W(λ⃗)
        λchain = √(I₁(λ⃗) / 3)
        W8(sum(λ⃗) / √(3N) - λchain / √(N)) + W8(λchain / √(N))
    end
end

"""
Davidson - Goulbourne

Parameters: Gc, Ge, λmax

Model: ``\\frac{G_c}{6}I_1-G_c\\lambda_{max}\\log\\bigg(3\\lambda_{max}^2-I_1\\bigg)+G_e\\sum\\limits_{i=1}^{3}\\big(\\lambda_i+\\frac{1}{\\lambda_i}\\big)``
"""
function DavidsonGoulbourne((; Gc, Ge, λmax))
    W(λ⃗) = 1 / 6 * Gc * I₁(λ⃗) - Gc * λmax^2 * log(3λmax^2 - I₁(λ⃗)) + Ge * sum(λ⃗ .+ 1 ./ λ⃗)
end

"""
Network Averaging Tube

Parameters: μcκ, n, q, μt

Model: ``\\mu_c \\kappa n \\log\\bigg(\\frac{\\sin(\\frac{\\pi}{\\sqrt{n}})(\\frac{I_1}{3})^{\\frac{q}{2}}}{\\sin(\\frac{\\pi}{\\sqrt{n}}(\\frac{I_1}{3})^{\\frac{q}{2}}}+\\mu_t\\big[\\frac{I_2}{3}^{1/2} - 1 \\big]``
"""
function NetworkAveragingTube((; μcκ, n, q, μt))
    W(λ⃗) = μcκ * n * log((sin(π / sqrt(n)) * (I₁(λ⃗) / 3)^(q / 2)) / (sin(π / sqrt(n) * (I₁(λ⃗) / 3)^(q / 2)))) + μt * ((I₂(λ⃗) / 3)^(1 / 2) - 1)
end

"""
General Constitutive Model

Parameters: Gc, Ge, N

Model: ``G_c N \\log\\bigg(\\frac{3N+\\frac{1}{2}I_1}{3N-I_1}\\bigg)+G_e\\sum\\limits_{i=1}^{3}\\frac{1}{\\lambda_I}``
"""
function GeneralConstitutiveModel((;Gc, Ge, N))
    W(λ⃗) = Gc*N*log((3N+0.5*I₁(λ⃗))/(3N-I₁(λ⃗))) + Ge*sum(λ⃗.^(-1))
end

"""
Full Network - Wu Geisson

Parameters: μ, N, ρ

Model: ``(1-\\rho)W_{3Chain}+\\rho W_{8chain}``
"""
function FullNetwork((;μ, N, ρ))
    W3 = ThreeChainModel((μ=μ, N=N))
    W8 = ArrudaBoyce((μ=μ, N=N))
    W(λ⃗) = (1-ρ)*W3(λ⃗)+ρ*W8(λ⃗)
end

"""
Zuniga - Beatty

Parameters: μ, N₃, N₈

Model: ``\\sqrt{\\frac{N_3+N_8}{2N_3}}W_{3Chain}+\\sqrt{\\frac{I_1}{3N_8}}W_{8Chain}``
"""
function ZunigaBeatty((;μ, N₃, N₈))
    ΛL  = √((N₃+N₈)/2)
    Λch = 1/√(3)*√(I₁(λ⃗))
    ρ₃  = ΛL/√(N₃)
    ρ₈  = Λch/√(N₈)
    W3 = ThreeChainModel((μ=μ, N=N₃))
    W8 = ArrudaBoyce((μ=μ, N₈=N₈))
    W(λ⃗) = ρ₃ * W3(λ⃗) + ρ₈ * W8(λ⃗)
end

"""
Lim

Parameters: μ₁, μ₂, N, Î₁

Model: ``(1-f(\\frac{I_1-3}{\\hat{I_1}-3}))W_{NeoHookean}(μ₁)+fW_{ArrudaBoyce}(μ₂, N)``
"""
function Lim((;μ₁, μ₂, N, Î₁); f = (x)->x^3*(10-15x+6x^2))
    Wg = NeoHookean((μ=μ₁))
    W8 = ArrudaBoyce((μ=μ₂, N = N))
    function W(λ⃗)
        ζ = (I₁-3)/(Î₁-3)
        (1-f(ζ))*Wg(λ⃗)+f(ζ)*W8(λ⃗)
    end
end

"""
Bechir Chevalier

Parameters: μ₀, η, ρ, N₃, N₈

Model: 

``W_{3Chain}(\\mu_f, N_3)+W_{8Chain}(\\frac{\\mu_c}{3}, N_8)``

``\\mu_f = \\rho\\sqrt{\\frac{I_1}{3N_8}}``

``\\mu_c = \\bigg(1-\\frac{\\eta\\alpha}{\\sqrt{N_3}}\\bigg)\\mu_0``

``\\alpha = \\max{\\lambda_1, \\lambda_2, \\lambda_3}``
"""
function BechirChevalier((;μ₀, η, ρ, N₃, N₈))
    μf = ρ*√(I₁/3/N₈)
    W3 = ThreeChainModel((μ = μf,N= N₃))
    function W(λ⃗)
        α = maximum(λ⃗)
        μc = (1-η*α/√(N₃))*μ₀
        W8 = ArrudaBoyce((μ = μc/3, N=N₈))
        W3(λ⃗)+W8(λ⃗)
    end
end