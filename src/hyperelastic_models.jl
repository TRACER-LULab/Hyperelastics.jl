export GeneralizedMooneyRivlin, GeneralizedDarijaniNaghdabadi, GeneralizedBeda

export MooneyRivlin, NeoHookean, Gent, Biderman, Isihara, JamesGreenSimpson, Lion, Yeoh, HauptSedlan, HartmannNeff, HainesWilson, Carroll, Nunes, BahremanDarijani, Zhao

export Knowles, Swanson, YamashitaKawabata, DavisDeThomas, Gregory, ModifiedGregory, Beda, Amin, LopezPamies, GenYeoh, HartSmith, VerondaWestmann, FungDemiray, Vito, HumphreyYin, ModifiedYeoh, Martins, ChevalierMarco, GornetDesmorat, MansouriDarijani, GentThomas, Alexander, LambertDianiRey, HossMarczakI, HossMarczakII, ExpLn

export Kilian, VanDeWaals, TakamizawaHayashi, YeohFleming, Gent3, PucciSaccomandi, HorganSaccomandi, Beatty, HorganMurphy

## Table 1
function GeneralizedMooneyRivlin((; C))
    function (λ⃗)
        I1 = I₁(λ⃗)
        I2 = I₂(λ⃗)
        @tullio W := C[j, i] * (I1 - 3)^(i - 1) * (I2 - 3)^(j - 1)
        return W
    end
end

function GeneralizedDarijaniNaghdabadi((;A, B, m, n))
    (λ⃗) -> sum(A .* (λ⃗ .^ m .- 1) + B .* (λ⃗ .^ (-n) .- 1))
end

function GeneralizedBeda((;C, K, α, β))
    function (λ⃗)
        # @tullio W:= C[i]/α[i]*(I₁(λ⃗)-3)^α[i]+K[j]/β[j]*(I₂(λ⃗)-3)^β[j]
        W1 = C ./ α .* (I₁(λ⃗) - 3) .^ α  |> sum
        W2 = K ./ β .* (I₂(λ⃗) - 3) .^ β |> sum
        return W1+W2
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
            A = [0.0, A2, 0.0, A4, 0.0, A6], 
            B = [0.0, B2], 
            m = [0.0, 2.0, 0.0, 4.0, 0.0, 6.0], 
            n = [0.0, 2.0])
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
        C = [C1, C2, C3],
        K = [K1],
        α = [α, 1.0, ζ],
        β = [β]
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


function GenYeoh((; K1, K2, K3, m, p, q))
    (λ⃗) -> K1 * (I₁(λ⃗) - 3)^m + K2 * (I₁(λ⃗) - 3)^p + K3 * (I₁(λ⃗) - 3)^q
end

function VerondaWestmann((; C1, C2, α))
    (λ⃗) -> C1 * (exp(α * (I₁(λ⃗) - 3)) - 1) + C2 * (I₂(λ⃗) - 3)
end

function FungDemiray((; μ, b))
    (λ⃗) -> μ / (2 * b) * (exp(b * (I₁(λ⃗) - 3)) - 1)
end

function Vito((; μ, b, α))
    (λ⃗) -> μ / (2 * b) * (exp(b * (α * (I₁(λ⃗) - 3) + (1 - α) * (I₂(λ⃗) - 3))) - 1)
end

function HumphreyYin((; C1, C2))
    (λ⃗) -> C1 * (exp(C2 * (I₁(λ⃗) - 3)) - 1)
end

function ModifiedYeoh((; C10, C20, C30, α, β))
    (λ⃗) -> C10 * (I₁(λ⃗) - 3) + C20 * (I₁(λ⃗) - 3)^2 + C30 * (I₁(λ⃗) - 3)^3 + α / β * (1 - exp(-β * (I₁(λ⃗) - 3)))
end

function MansouriDarijani((; A1, m1, B1, n1))
    (λ⃗) -> A1 * (exp(m1 * (I₁(λ⃗) - 3)) - 1) + B1 * (exp(n1 * (I₂(λ⃗) - 3)) - 1)
end

function GentThomas((; C1, C2))
    (λ⃗) -> C1 * (I₁(λ⃗) - 3) + C2 * log(I₂(λ⃗) / 3)
end

function HossMarczakI((; α, β, μ, b, n))
    (λ⃗) -> α / β * (1 - exp(-β * (I₁(λ⃗) - 3))) + μ / (2b) * ((1 + b / n * (I₁(λ⃗) - 3))^n - 1)
end

function HossMarczakII((; α, β, μ, b, n, C2))
    (λ⃗) -> α / β * (1 - exp(-β * (I₁(λ⃗) - 3))) + μ / (2b) * ((1 + b / n * (I₁(λ⃗) - 3))^n - 1) + C2 * log(I₂(λ⃗) / 3)
end

function ExpLn((; A, a, b))
    (λ⃗) -> A * (1 / a * exp(a * (I₁(λ⃗) - 3)) + b * (I₁(λ⃗) - 2) * (1 - log(I₁(λ⃗) - 2)) - 1 / a - b)
end

function Killian((; μ, JL))
    (λ⃗) -> -μ * JL * (log(1 - sqrt((I₁(λ⃗) - 3) / JL)) + sqrt((I₁(λ⃗) - 3) / JL))
end

# function VanDerWaals((;μ, λm, β, α))
#     (λ⃗) -> μ*(-(λm^2-3)*(log(1-θ)+θ)-2/3*α*((I₁(λ⃗)-3)/2)^(3/2))
# end

# function TakamizawaHayashi((;c, Jm))
#     (λ⃗) -> -c*log(1-((I₁(λ⃗)-3)/Jm)^2)
# end

function YeohFleming((;A, B, C10, Im))
    (λ⃗) -> A/B*(1-exp(-B*(I₁(λ⃗)-3)))-C10*(Im-3)*log(1-((I₁(λ⃗)-3)/(Im-3)))
end


# ---------------------------------------------------- #
function Gent((; μ, Jₘ))
    (λ⃗) -> -μ * Jₘ / 2 * log(1 - (I₁(λ⃗) - 3) / Jₘ)
end
