# # Available Models
export GeneralMooneyRivlin, GeneralDarijaniNaghdabadi, GeneralBeda, MooneyRivlin, NeoHookean, Gent, Biderman, Isihara, JamesGreenSimpson, Lion, Yeoh, HauptSedlan, HartmannNeff, HainesWilson, Carroll, BahremanDarijani, Zhao, Knowles, Swanson, YamashitaKawabata, DavisDeThomas, Gregory, ModifiedGregory, Beda, Amin, LopezPamies, GenYeoh, VerondaWestmann, FungDemiray, Vito, ModifiedYeoh, MansouriDarijani, GentThomas, HossMarczakI, HossMarczakII, ExpLn, VanDerWaals, TakamizawaHayashi, YeohFleming, PucciSaccomandi, HorganSaccomandi, Beatty, HorganMurphy, ArrudaBoyce, Ogden, EdwardVilgis, NonaffineTube, Tube, MCC, Bechir4Term, ConstrainedJunction, ContinuumHybrid, ArmanNarooei, PengLandel, ValanisLandel, Attard, Shariff, ThreeChainModel, ModifiedFloryErman, ABGI, BechirChevalier, Bootstrapped8Chain, DavidsonGoulbourne, ExtendedTubeModel, FullNetwork, HartSmith, GeneralConstitutiveModel, Lim, NonaffineMicroSphere, AffineMicroSphere, KhiemItskov, ZunigaBeatty, ChevalierMarco, Alexander, GornetDesmorat, LambertDianiRey

"""
General Mooney Rivlin[^1]

Parameters: [C]

Model:
``\\sum\\limits_{i,j = 0}^{N,M} C_{i,j}(I_1-3)^i(2I_2-3)^j``

[^1]: > Mooney M. A theory of large elastic deformation. Journal of applied physics. 1940 Sep;11(9):582-92.
"""
struct GeneralMooneyRivlin end

function StrainEnergyDensityFunction(ψ::GeneralMooneyRivlin, (; C))
    function W(λ⃗)
        I1 = I₁(λ⃗)
        I2 = I₂(λ⃗)
        @tullio W := C[j, i] * (I1 - 3)^(i - 1) * (I2 - 3)^(j - 1)
        return W
    end
end

function StrainEnergyDensityFunction(ψ::GeneralMooneyRivlin, (; C), I::InvariantForm)
    function W(I⃗)
        I₁, I₂ = I⃗
        @tullio W := C[j, i] * (I₁ - 3)^(i - 1) * (I₂ - 3)^(j - 1)
        return W
    end
end

function parameters(ψ::GeneralMooneyRivlin)
    return (:μ, :Jₘ)
end

"""
General Darijani Naghdabadi [^1]

Parameters: A⃗, B⃗, m⃗, n⃗

Model: ``\\sum\\limits_{i = 1}{3}\\sum\\limits_{j=0}^{N} A_j (\\lambda_i^{m_j}-1) + B_j(\\lambda_i^{-n_j}-1)``

[^1]: > Bahreman M, Darijani H. New polynomial strain energy function; application to rubbery circular cylinders under finite extension and torsion. Journal of Applied Polymer Science. 2015 Apr 5;132(13).
"""
struct GeneralDarijaniNaghdabadi end


function StrainEnergyDensityFunction(ψ::GeneralDarijaniNaghdabadi, (; A⃗, B⃗, m⃗, n⃗))
    @assert length(A⃗)==length(B⃗)==length(m⃗)==length(n⃗) "The vectors are not the same length"
    W(λ⃗) = sum(A⃗ .* (λ⃗ .^ m⃗ .- 1) + B⃗ .* (λ⃗ .^ (-n⃗) .- 1))
end

function parameters(ψ::GeneralDarijaniNaghdabadi)
    return (:A⃗, :B⃗, :m⃗, :n⃗)
end


"""
General Beda [^1]

Parameters: C, K, α, β

Model: ``\\sum\\limits_{i = 1}^{N}\\frac{C_i}{\\alpha_i}(I_1-3)^{\\alpha_i} + \\sum\\limits_{j=1}^{M}\\frac{K_j}{\\beta_j}(I_2-3)^{\\beta_j}``

[^1]: > Beda T. Reconciling the funda
mental phenomenological expression of the strain energy of rubber with established experimental facts. Journal of Polymer Science Part B: Polymer Physics. 2005 Jan 15;43(2):125-34.
"""
struct GeneralBeda end

function StrainEnergyDensityFunction(ψ::GeneralBeda, (; C, K, α, β))
    @assert length(C)==length(α) "Vector C and Vector α are not the same length"
    @assert length(K)==length(β) "Vector K and Vector β are not the same length"
    function W(λ⃗)
        W1 = C ./ α .* (I₁(λ⃗) - 3) .^ α |> sum
        W2 = K ./ β .* (I₂(λ⃗) - 3) .^ β |> sum
        return W1 + W2
    end
end

function StrainEnergyDensityFunction(ψ::GeneralBeda, (; C, K, α, β), I::InvariantForm)
    @assert length(C)==length(α) "Vector C and Vector α are not the same length"
    @assert length(K)==length(β) "Vector K and Vector β are not the same length"
    function W(I⃗)
        W1 = C ./ α .* (I⃗[1] - 3) .^ α |> sum
        W2 = K ./ β .* (I⃗[2] - 3) .^ β |> sum
        return W1 + W2
    end
end

function parameters(ψ::GeneralBeda)
    return (:C, :K, :α, :β)
end

"""
Mooney Rivlin Model [^1]

Parameters: C01, C10

Model: ``C_{10}(I_1-3)+C_{01}(I_2-3)``

[^1]: > Mooney M. A theory of large elastic deformation. Journal of applied physics. 1940 Sep;11(9):582-92.
"""
struct MooneyRivlin end

function StrainEnergyDensityFunction(ψ::MooneyRivlin, (; C10, C01))
    StrainEnergyDensityFunction(
        GeneralMooneyRivlin(),
        (C=[
            0.0 C10
            C01 0.0
        ],
        )
    )
end

function StrainEnergyDensityFunction(ψ::MooneyRivlin, (; C10, C01), I::InvariantForm)
    StrainEnergyDensityFunction(
        GeneralMooneyRivlin(),
        (C=[
            0.0 C10
            C01 0.0
        ],
        ),
        I
    )
end

function parameters(ψ::MooneyRivlin)
    return (:C10, :C01)
end


"""
NeoHookean [^1]

Parameters: μ

Model: ``\\frac{\\mu}{2}(I_1-3)``

[^1]: > Treloar LR. The elasticity of a network of long-chain molecules—II. Transactions of the Faraday Society. 1943;39:241-6.
"""
struct NeoHookean end

function StrainEnergyDensityFunction(ψ::NeoHookean, (; μ))
    W(λ⃗) = μ / 2 * (I₁(λ⃗) - 3)
end

function StrainEnergyDensityFunction(ψ::NeoHookean, (; μ), I::InvariantForm)
    function W(I⃗)
        μ / 2 * (I⃗[1] - 3)
    end
end

function parameters(ψ::NeoHookean)
    return (:μ)
end

"""
Isihara [^1]

Parameters: C10, C20, C01

Model: ``\\sum\\limits_{i,j=0}^{2, 1}C_{i,j}(I_1-3)^i(I_2-3)^j``

[^1]: > Isihara A, Hashitsume N, Tatibana M. Statistical theory of rubber‐like elasticity. IV.(two‐dimensional stretching). The Journal of Chemical Physics. 1951 Dec;19(12):1508-12.
"""
struct Isihara end

function StrainEnergyDensityFunction(ψ::Isihara, (; C10, C20, C01))
    StrainEnergyDensityFunction(
        GeneralMooneyRivlin(),
        (C=[
            0.0 C10 C20
            C01 0.0 0.0
        ],
        )
    )
end

function StrainEnergyDensityFunction(ψ::Isihara, (; C10, C20, C01), I::InvariantForm)
    StrainEnergyDensityFunction(
        GeneralMooneyRivlin(),
        (C=[
            0.0 C10 C20
            C01 0.0 0.0
        ],
        ),
        I
    )
end

function parameters(ψ::Isihara)
    return (:C10, :C20, :C01)
end

"""
Biderman [^1]

Parameters: C10, C01, C20, C30

Model: ``\\sum\\limits_{i,j=0}^{3, 1}C_{i,j}(I_1-3)^i(I_2-3)^j``

[^1]: > Biderman VL. Calculation of rubber parts. Rascheti na prochnost. 1958;40.
"""
struct Biderman end

function StrainEnergyDensityFunction(ψ::Biderman, (; C10, C01, C20, C30))
    StrainEnergyDensityFunction(
        GeneralMooneyRivlin(),
        (C=[
            0.0 C10 C20 C30
            C01 0.0 0.0 0.0
        ],
        )
    )
end

function StrainEnergyDensityFunction(ψ::Biderman, (; C10, C01, C20, C30), I::InvariantForm)
    StrainEnergyDensityFunction(
        GeneralMooneyRivlin(),
        (C=[
            0.0 C10 C20 C30
            C01 0.0 0.0 0.0
        ],
        ),
        I
    )
end

function parameters(ψ::Biderman)
    return (:C10, :C01, :C20, :C30)
end

"""
James-Green-Simpson [^1]

Parameters: C10, C01, C11, C20, C30

Model: ``\\sum\\limits_{i,j=0}^{3, 1}C_{i,j}(I_1-3)^i(I_2-3)^j``

[^1]: > James AG, Green A, Simpson GM. Strain energy functions of rubber. I. Characterization of gum vulcanizates. Journal of Applied Polymer Science. 1975 Jul;19(7):2033-58.
"""
struct JamesGreenSimpson end

function StrainEnergyDensityFunction(ψ::JamesGreenSimpson, (; C10, C01, C11, C20, C30))
    StrainEnergyDensityFunction(
        GeneralMooneyRivlin(),
        (C=[
            0.0 C10 C20 C30
            C01 0.0 0.0 0.0
        ],
        )
    )
end

function StrainEnergyDensityFunction(ψ::JamesGreenSimpson, (; C10, C01, C11, C20, C30), I::InvariantForm)
    StrainEnergyDensityFunction(
        GeneralMooneyRivlin(),
        (C=[
            0.0 C10 C20 C30
            C01 0.0 0.0 0.0
        ],
        ),
        I
    )
end

function parameters(ψ::JamesGreenSimpson)
    return (:C10, :C01, :C11, :C20, :C30)
end

"""
Haines-Wilson [^1]

Parameters: C10, C01, C11, C02, C20, C30

Model: ``\\sum\\limits_{i,j=0}^{3, 2}C_{i,j}(I_1-3)^i(I_2-3)^j``

[^1]: > Haines DW, Wilson WD. Strain-energy density function for rubberlike materials. Journal of the Mechanics and Physics of Solids. 1979 Aug 1;27(4):345-60.
"""
struct HainesWilson end

function StrainEnergyDensityFunction(ψ::HainesWilson, (; C10, C01, C11, C02, C20, C30))
    StrainEnergyDensityFunction(
        GeneralMooneyRivlin(),
        (C=[
            0.0 C10 C20 C30
            C01 C11 0.0 0.0
            C02 0.0 0.0 0.0
        ],
        )
    )
end

function StrainEnergyDensityFunction(ψ::HainesWilson, (; C10, C01, C11, C02, C20, C30), I::InvariantForm)
    StrainEnergyDensityFunction(
        GeneralMooneyRivlin(),
        (C=[
            0.0 C10 C20 C30
            C01 C11 0.0 0.0
            C02 0.0 0.0 0.0
        ],
        ),
        I
    )
end

function parameters(ψ::HainesWilson)
    return (:C10, :C01, :C11, :C02, :C20, :C30)
end

"""
Yeoh [^1]

Parameters: C10, C20, C30

Model: ``\\sum\\limits_{i,j=0}^{3, 0}C_{i,j}(I_1-3)^i(I_2-3)^j``

[^1]: > Yeoh OH. Characterization of elastic properties of carbon-black-filled rubber vulcanizates. Rubber chemistry and technology. 1990 Nov;63(5):792-805.
"""
struct Yeoh end

function StrainEnergyDensityFunction(ψ::Yeoh, (; C10, C20, C30))
    StrainEnergyDensityFunction(
        GeneralMooneyRivlin(),
        (C=[0.0 C10 C20 C30],)
    )
end

function StrainEnergyDensityFunction(ψ::Yeoh, (; C10, C20, C30), I::InvariantForm)
    StrainEnergyDensityFunction(
        GeneralMooneyRivlin(),
        (C=[0.0 C10 C20 C30],),
        I
    )
end

function parameters(ψ::Yeoh)
    return (:C10, :C20, :C30)
end

"""
Lion [^1]

Parameters: C10, C01, C50

Model: ``\\sum\\limits_{i,j=0}^{5,1}C_{i,j}(I_1-3)^i(I_2-3)^j``

[^1]: > Lion A. On the large deformation behaviour of reinforced rubber at different temperatures. Journal of the Mechanics and Physics of Solids. 1997 Nov 1;45(11-12):1805-34.
"""
struct Lion end

function StrainEnergyDensityFunction(ψ::Lion, (; C10, C01, C50))
    StrainEnergyDensityFunction(
        GeneralMooneyRivlin(),
        (C=[
            0.0 C10 0.0 0.0 0.0 C50
            C01 0.0 0.0 0.0 0.0 0.0
        ],)
    )
end

function StrainEnergyDensityFunction(ψ::Lion, (; C10, C01, C50), I::InvariantForm)
    StrainEnergyDensityFunction(
        GeneralMooneyRivlin(),
        (C=[
            0.0 C10 0.0 0.0 0.0 C50
            C01 0.0 0.0 0.0 0.0 0.0
        ],),
        I
    )
end

function parameters(ψ::Lion)
    return (:C10, :C01, :C50)
end


"""
Haupt Sedlan [^1]

Parameters: C10, C01, C11, C02, C30

Model:
``\\sum\\limits_{i,j=0}^{3, 2}C_{i,j}(I_1-3)^i(I_2-3)^j``

[^1]: > Haupt P, Sedlan K. Viscoplasticity of elastomeric materials: experimental facts and constitutive modelling. Archive of Applied Mechanics. 2001 Mar;71(2):89-109.
"""
struct HauptSedlan end

function StrainEnergyDensityFunction(ψ::HauptSedlan, (; C10, C01, C11, C02, C30))
    StrainEnergyDensityFunction(
        GeneralMooneyRivlin(),
        (C=[
            0.0 C10 0.0 C30
            C01 C11 0.0 0.0
            C02 0.0 0.0 0.0
        ],)
    )
end

function StrainEnergyDensityFunction(ψ::HauptSedlan, (; C10, C01, C11, C02, C30), I::InvariantForm)
    StrainEnergyDensityFunction(
        GeneralMooneyRivlin(),
        (C=[
            0.0 C10 0.0 C30
            C01 C11 0.0 0.0
            C02 0.0 0.0 0.0
        ],),
        I
    )
end

function parameters(ψ::HauptSedlan)
    return (:C10, :C01, :C11, :C02, :C30)
end

"""
Hartmann-Neff [^1]

Parameters: α, Ci0, C0j

Model: ``\\sum\\limits_{i,j=0}^{M,N}C_{i,0}(I_1-3)^i -3\\sqrt{3}^j+\\alpha(I_1-3)``

[^1]: > Hartmann S, Neff P. Polyconvexity of generalized polynomial-type hyperelastic strain energy functions for near-incompressibility. International journal of solids and structures. 2003 Jun 1;40(11):2767-91.
"""
struct HartmannNeff end

function StrainEnergyDensityFunction(ψ::HartmannNeff, (; α, Ci0, C0j))
    function W(λ⃗)
        @tullio W1 := Ci0[i] * (I₁(λ⃗) - 3)^i
        @tullio W2 := C0j[j] * (I₂(λ⃗)^(3 / 2) - 3sqrt(3))^j
        return W1 + W2 + α * (I₁(λ⃗)^3 - 3^3)
    end
end

function StrainEnergyDensityFunction(ψ::HartmannNeff, (; α, Ci0, C0j), I::InvariantForm)
    function W(I⃗)
        I₁, I₂ = I⃗
        @tullio W1 := Ci0[i] * (I₁ - 3)^i
        @tullio W2 := C0j[j] * (I₂^(3 / 2) - 3sqrt(3))^j
        return W1 + W2 + α * (I₁^3 - 3^3)
    end
end

function parameters(ψ::HartmannNeff)
    return (:α, :Ci0, :C0j)
end

"""
Carroll [^1]

Parameters: A, B, C

Model: ``AI_1+BI_1^4+C\\sqrt{I_2}``

[^1]: > Carroll M. A strain energy function for vulcanized rubbers. Journal of Elasticity. 2011 Apr;103(2):173-87.
"""
struct Carroll end

function StrainEnergyDensityFunction(ψ::Carroll, (; A, B, C))
    W(λ⃗) = A * I₁(λ⃗) + B * I₁(λ⃗)^4 + C * I₂(λ⃗)^(1 / 2)
end

function StrainEnergyDensityFunction(ψ::Carroll, (; A, B, C), I::InvariantForm)
    W(I⃗) = A * I⃗[1] + B * I⃗[1]^4 + C * I⃗[2]^(1 / 2)
end

function parameters(ψ::Carroll)
    return (:A, :B, :C)
end

"""
Bahreman Darijani [^1]

Parameters: A2, B2, A4, A6

Model:
``\\sum\\limits_{i = 1}{3}\\sum\\limits_{j=0}^{N} A_j (\\lambda_i^{m_j}-1) + B_j(\\lambda_i^{-n_j}-1)``

[^1]: > Bahreman M, Darijani H. New polynomial strain energy function; application to rubbery circular cylinders under finite extension and torsion. Journal of Applied Polymer Science. 2015 Apr 5;132(13).
"""
struct BahremanDarijani end

function StrainEnergyDensityFunction(ψ::BahremanDarijani, (; A2, B2, A4, A6))
    StrainEnergyDensityFunction(
        GeneralDarijaniNaghdabadi(),
        (
            A=[0.0, A2, 0.0, A4, 0.0, A6],
            B=[0.0, B2],
            m=[0.0, 2.0, 0.0, 4.0, 0.0, 6.0],
            n=[0.0, 2.0])
    )
end

function parameters(ψ::BahremanDarijani)
    return (:A2, :B2, :A4, :A6)
end

"""
Zhao [^1]

Parameters: C₋₁¹,, C₁¹, C₂¹, C₂²

Model: ``C_{-1}^1*(I_2-3)+C_{1}^{1}(I_1-3)+C_{2}^{1}(I_1^2-2I_2-3)+C_{2}^{2}(I_1^2-2I_2-3)^2``

[^1]: > Zhao Z, Mu X, Du F. Modeling and verification of a new hyperelastic model for rubber-like materials. Mathematical Problems in Engineering. 2019 May 2;2019.
"""
struct Zhao end

function StrainEnergyDensityFunction(ψ::Zhao, (; C₋₁¹, C₁¹, C₂¹, C₂²))
    W(λ⃗) = C₋₁¹ * (I₂(λ⃗) - 3) + C₁¹ * (I₁(λ⃗) - 3) + C₂¹ * (I₁(λ⃗)^2 - 2I₂(λ⃗) - 3) + C₂² * (I₁(λ⃗)^2 - 2I₂(λ⃗) - 3)^2
end

function StrainEnergyDensityFunction(ψ::Zhao, (; C₋₁¹, C₁¹, C₂¹, C₂²))
    W(I⃗) = C₋₁¹ * (I⃗[2] - 3) + C₁¹ * (I⃗[1] - 3) + C₂¹ * (I⃗[1]^2 - 2I⃗[2] - 3) + C₂² * (I⃗[1]^2 - 2I⃗[2] - 3)^2
end

function parameters(ψ::Zhao)
    return (:C₋₁¹, :C₁¹, :C₂¹, :C₂²)
end

"""
Knowles [^1]

Parameters: μ, b, n

Model: ``\\frac{\\mu}{2b}((1+\\frac{b}{n}(I_1-3))^n-1)``

[^1]: > Knowles JK. The finite anti-plane shear field near the tip of a crack for a class of incompressible elastic solids. International Journal of Fracture. 1977 Oct;13(5):611-39.
"""
struct Knowles end

function StrainEnergyDensityFunction(ψ::Knowles, (; μ, b, n))
    W(λ⃗) = μ / (2b) * ((1 + (b / n) * (I₁(λ⃗) - 3))^n - 1)
end

function StrainEnergyDensityFunction(ψ::Knowles, (; μ, b, n), I::InvariantForm)
    W(I⃗) = μ / (2b) * ((1 + (b / n) * (I⃗[1] - 3))^n - 1)
end


function parameters(ψ::Knowles)
    return (:μ, :b, :n)
end

"""
Swanson [^1]

Parameters: A, α, B, β

Model: ``\\sum\\limits_{i=1}^{N} \\frac{3}{2}(\\frac{A_i}{1+\\alpha_i}(\\frac{I_1}{3})^{1+\\alpha_i}+\\frac{B_i}{1+\\beta_i}(\\frac{I_2}{3})^{1+\\beta_i}``

[^1]: > Swanson SR. A constitutive model for high elongation elastic materials.
"""
struct Swanson end

function StrainEnergyDensityFunction(ψ::Swanson, (; A, α, B, β))
    @assert length(A) == length(α) == length(B) == length(β) "The vectors are not the same length"
    W(λ⃗) = @tullio _ := 3 / 2 * (A[i] / (1 + α[i]) * (I₁(λ⃗) / 3)^(1 + α[i]) + B[i] / (1 + β[i]) * (I₂(λ⃗) / 3)^(1 + β[i]))
end

function StrainEnergyDensityFunction(ψ::Swanson, (; A, α, B, β), I::InvariantForm)
    @assert length(A) == length(α) == length(B) == length(β) "The vectors are not the same length"
    W(I⃗) = @tullio _ := 3 / 2 * (A[i] / (1 + α[i]) * (I⃗[1] / 3)^(1 + α[i]) + B[i] / (1 + β[i]) * (I⃗[2] / 3)^(1 + β[i]))
end


function parameters(ψ::Swanson)
    return (:A, :α, :B, :β)
end

"""
Yamashita-Kawabata [^1]

Parameters: C1, C2, C3, N

Model: ``C_1(I_1-3)+C_2(I_2-3)+\\frac{C_3}{N+1}(I_1-3)^{N+1}``

[^1]: > Yamashita Y, Kawabata S. Approximated form of the strain energy-density function of carbon-black filled rubbers for industrial applications. Nippon Gomu Kyokaishi(Journal of the Society of Rubber Industry, Japan)(Japan). 1992;65(9):517-28.
"""
struct YamashitaKawabata end

function StrainEnergyDensityFunction(ψ::YamashitaKawabata, (; C1, C2, C3, N))
    W(λ⃗) = C1 * (I₁(λ⃗) - 3) + C2 * (I₂(λ⃗) - 3) + C3 / (N + 1) * (I₁(λ⃗) - 3)^(N + 1)
end

function StrainEnergyDensityFunction(ψ::YamashitaKawabata, (; C1, C2, C3, N), I::InvariantForm)
    W(I⃗) = C1 * (I⃗[1] - 3) + C2 * (I⃗[2] - 3) + C3 / (N + 1) * (I⃗[1] - 3)^(N + 1)
end

function parameters(ψ::YamashitaKawabata)
    return (:C1, :C2, :C3, :N)
end

"""
Davis-DeThomas [^1]

Parameters: A, n, C, k

Model: ``\\frac{A}{2(1-\\frac{n}{2})}(I_1-3+C^2)^{1-\\frac{n}{2}}+k(I_1-3)^2``

[^1]: > Davies CK, De DK, Thomas AG. Characterization of the behavior of rubber for engineering design purposes. 1. Stress-strain relations. Rubber chemistry and technology. 1994 Sep;67(4):716-28.
"""
struct DavisDeThomas end

function StrainEnergyDensityFunction(ψ::DavisDeThomas, (; A, n, C, k))
    W(λ⃗) = A / (2 * (1 - n / 2)) * (I₁(λ⃗) - 3 + C^2)^(1 - n / 2) + k * (I₁(λ⃗) - 3)^2
end

function StrainEnergyDensityFunction(ψ::DavisDeThomas, (; A, n, C, k))
    W(I⃗) = A / (2 * (1 - n / 2)) * (I⃗[1] - 3 + C^2)^(1 - n / 2) + k * (I⃗[1] - 3)^2
end

function parameters(ψ::DavisDeThomas)
    return (:A, :n, :C, :k)
end

"""
Gregory [^1]

Parameters: A, B, C, m, n

Model: ``\\frac{A}{2-n}(I_1-3+C^2)^{1-\\frac{n}{2}}+\\frac{B}{2+m}(I_1-3+C^2)^{1+\\frac{m}{2}}``

[^1]: > Gregory IH, Muhr AH, Stephens IJ. Engineering applications of rubber in simple extension. Plastics rubber and composites processing and applications. 1997;26(3):118-22.
"""
struct Gregory end

function StrainEnergyDensityFunction(ψ::Gregory, (; A, B, C, m, n))
    W(λ⃗) = A / (2 - n) * (I₁(λ⃗) - 3 + C^2)^(1 - n / 2) + B / (2 + m) * (I₁(λ⃗) - 3 + C^2)^(1 + m / 2)
end

function StrainEnergyDensityFunction(ψ::Gregory, (; A, B, C, m, n), I::InvariantForm)
    W(I⃗) = A / (2 - n) * (I⃗[1] - 3 + C^2)^(1 - n / 2) + B / (2 + m) * (I⃗[1] - 3 + C^2)^(1 + m / 2)
end

function parameters(ψ::Gregory)
    return (:A, :B, :C, :m, :n)
end

"""
Modified Gregory [^1]

Parameters: A, α, M, B, β, N

Model: ``\\frac{A}{1+\\alpha}(I_1-3+M^2)^{1+\\alpha}+\\frac{B}{1+\\beta}(I_1-3+N^2)^{1+\\beta}``

[^1]: > He H, Zhang Q, Zhang Y, Chen J, Zhang L, Li F. A comparative study of 85 hyperelastic constitutive models for both unfilled rubber and highly filled rubber nanocomposite material. Nano Materials Science. 2021 Jul 16.
"""
struct ModifiedGregory end

function StrainEnergyDensityFunction(ψ::ModifiedGregory, (; A, α, M, B, β, N))
    W(λ⃗) = A / (1 + α) * (I₁(λ⃗) - 3 + M^2)^(1 + α) + B / (1 + β) * (I₁(λ⃗) - 3 + N^2)^(1 + β)
end

function StrainEnergyDensityFunction(ψ::ModifiedGregory, (; A, α, M, B, β, N))
    W(I⃗) = A / (1 + α) * (I⃗[1] - 3 + M^2)^(1 + α) + B / (1 + β) * (I⃗[1] - 3 + N^2)^(1 + β)
end

function parameters(ψ::ModifiedGregory)
    return (:A, :α, :M, :B, :β, :N)
end

"""
Beda [^1]

Parameters: C1, C2, C3, K1, α, β, ζ

Model: ``\\frac{C_1}{\\alpha}(I_1-3)^{\\alpha}+C_2(I_1-3)+\\frac{C_3}{\\zeta}(I_1-3)^{\\zeta}+\\frac{K_1}{\\beta}(I_2-3)^\\beta``

[^1]: > Beda T. Reconciling the fundamental phenomenological expression of the strain energy of rubber with established experimental facts. Journal of Polymer Science Part B: Polymer Physics. 2005 Jan 15;43(2):125-34.
"""
struct Beda end

function StrainEnergyDensityFunction(ψ::Beda, (; C1, C2, C3, K1, α, β, ζ))
    StrainEnergyDensityFunction(
        GeneralBeda(),
        (
            C=[C1, C2, C3],
            K=[K1],
            α=[α, 1.0, ζ],
            β=[β]
        )
    )
end

function StrainEnergyDensityFunction(ψ::Beda, (; C1, C2, C3, K1, α, β, ζ), I::InvariantForm)
    StrainEnergyDensityFunction(
        GeneralBeda(),
        (
            C=[C1, C2, C3],
            K=[K1],
            α=[α, 1.0, ζ],
            β=[β]
        ),
        I
    )
end

function parameters(ψ::Beda)
    return (:C1, :C2, :C3, :K1, :α, :β, :ζ)
end

"""
Amin [^1]

Parameters: C1, C2, C3, C4, N, M

Model:``C_1 (I_1 - 3) + \\frac{C_2}{N + 1} (I_1 - 3)^{N + 1} + \\frac{C_3}{M + 1} (I_1 - 3)^{M + 1} + C_4 (I_2 - 3)``

[^1]: > Amin AF, Wiraguna SI, Bhuiyan AR, Okui Y. Hyperelasticity model for finite element analysis of natural and high damping rubbers in compression and shear. Journal of engineering mechanics. 2006 Jan;132(1):54-64.
"""
struct Amin end

function StrainEnergyDensityFunction(ψ::Amin, (; C1, C2, C3, C4, N, M))
    W(λ⃗) = C1 * (I₁(λ⃗) - 3) + C2 / (N + 1) * (I₁(λ⃗) - 3)^(N + 1) + C3 / (M + 1) * (I₁(λ⃗) - 3)^(M + 1) + C4 * (I₂(λ⃗) - 3)
end

function StrainEnergyDensityFunction(ψ::Amin, (; C1, C2, C3, C4, N, M), I::InvariantForm)
    W(I⃗) = C1 * (I⃗[1] - 3) + C2 / (N + 1) * (I⃗[1] - 3)^(N + 1) + C3 / (M + 1) * (I⃗[1] - 3)^(M + 1) + C4 * (I⃗[2] - 3)
end

function parameters(ψ::Amin)
    return (:C1, :C2, :C3, :C4, :N, :M)
end

"""
Lopez-Pamies [^1]

Parameters: α⃗, μ⃗

Model: ``\\frac{3.0^{1 - \\alpha_i}}{2\\alpha_i} \\mu_i (I_1^{\\alpha_i} - 3^{\\alpha_i})``

[^1]: > Lopez-Pamies O. A new I1-based hyperelastic model for rubber elastic materials. Comptes Rendus Mecanique. 2010 Jan 1;338(1):3-11.
"""
struct LopezPamies end

function StrainEnergyDensityFunction(ψ::LopezPamies, (; α⃗, μ⃗))
    @assert length(α⃗) == length(μ⃗) "length of α⃗ is not equal to length of μ⃗"
    W(λ⃗) = @tullio _ := (3.0^(1 - α⃗[i])) / (2α⃗[i]) * μ⃗[i] * (I₁(λ⃗)^(α⃗[i]) - 3^(α⃗[i]))
end

function StrainEnergyDensityFunction(ψ::LopezPamies, (; α⃗, μ⃗), I::InvariantForm)
    @assert length(α⃗) == length(μ⃗) "length of α⃗ is not equal to length of μ⃗"
    W(I⃗) = @tullio _ := (3.0^(1 - α⃗[i])) / (2α⃗[i]) * μ⃗[i] * (I⃗[1]^(α⃗[i]) - 3^(α⃗[i]))
end

function parameters(ψ::LopezPamies)
    return (:α⃗, :μ⃗)
end

"""
GenYeoh [^1]

Parameters: K1, K2, K3, m, p, q

Model: ``K_1 (I_1 - 3)^m + K_2 * (I_1 - 3)^p + K_3 * (I_1 - 3)^q``

[^1]: > Hohenberger TW, Windslow RJ, Pugno NM, Busfield JJ. A constitutive model for both low and high strain nonlinearities in highly filled elastomers and implementation with user-defined material subroutines in ABAQUS. Rubber Chemistry and Technology. 2019;92(4):653-86.
"""
struct GenYeoh end

function StrainEnergyDensityFunction(ψ::GenYeoh, (; K1, K2, K3, m, p, q))
    W(λ⃗) = K1 * (I₁(λ⃗) - 3)^m + K2 * (I₁(λ⃗) - 3)^p + K3 * (I₁(λ⃗) - 3)^q
end

function StrainEnergyDensityFunction(ψ::GenYeoh, (; K1, K2, K3, m, p, q))
    W(I⃗) = K1 * (I⃗[1] - 3)^m + K2 * (I⃗[1] - 3)^p + K3 * (I⃗[1] - 3)^q
end

function parameters(ψ::GenYeoh)
    return (:K1, :K2, :K3, :m, :p, :q)
end

"""
Hart-Smith [^1]

Parameters: G, k₁, k₂

Model: ``\\frac{G\\exp{(-9k_1+k_1I_1)}}{k_1}+Gk_2\\log{I_2}``

[^1]: > Hart-Smith LJ. Elasticity parameters for finite deformations of rubber-like materials. Zeitschrift für angewandte Mathematik und Physik ZAMP. 1966 Sep;17(5):608-26.
"""
struct HartSmith end

function StrainEnergyDensityFunction(ψ::HartSmith, (; G, k₁, k₂))
    W(λ⃗) = G * exp(-9k₁ + k₁ * I₁(λ⃗)) / k₁ + G * k₂ * log(I₂(λ⃗))
end

function StrainEnergyDensityFunction(ψ::HartSmith, (; G, k₁, k₂), I::InvariantForm)
    W(I⃗) = G * exp(-9k₁ + k₁ * I⃗[1]) / k₁ + G * k₂ * log(I⃗[2])
end

function parameters(ψ::HartSmith)
    return (:G, :k₁, :k₂)
end

"""
Veronda-Westmann [^1]

Parameters: C1, C2, α

Model: ``C_1 (\\exp(\\alpha(I_1 - 3)) - 1) + C_2 (I_2 - 3)``

[^1]: > Veronda DR, Westmann RA. Mechanical characterization of skin—finite deformations. Journal of biomechanics. 1970 Jan 1;3(1):111-24.
"""
struct VerondaWestmann end

function StrainEnergyDensityFunction(ψ::VerondaWestmann, (; C1, C2, α))
    W(λ⃗) = C1 * (exp(α * (I₁(λ⃗) - 3)) - 1) + C2 * (I₂(λ⃗) - 3)
end

function StrainEnergyDensityFunction(ψ::VerondaWestmann, (; C1, C2, α), I::InvariantForm)
    W(I⃗) = C1 * (exp(α * (I⃗[1] - 3)) - 1) + C2 * (I⃗[2] - 3)
end

function parameters(ψ::VerondaWestmann)
    return (:C1, :C2, :α)
end

"""
Fung-Demiray [^1][^2]

Parameters: μ, b

Model: ``\\frac{\\mu}{2 * b} (\\exp(b(I_1 - 3)) - 1)``

[^1]: > Fung YC. Elasticity of soft tissues in simple elongation. American Journal of Physiology-Legacy Content. 1967 Dec 1;213(6):1532-44.
[^2]: > Demiray H. A note on the elasticity of soft biological tissues. Journal of biomechanics. 1972 May 1;5(3):309-11.
"""
struct FungDemiray end

function StrainEnergyDensityFunction(ψ::FungDemiray, (; μ, b))
    W(λ⃗) = μ / (2 * b) * (exp(b * (I₁(λ⃗) - 3)) - 1)
end

function StrainEnergyDensityFunction(ψ::FungDemiray, (; μ, b), I::InvariantForm)
    W(I⃗) = μ / (2 * b) * (exp(b * (I⃗[1] - 3)) - 1)
end

function parameters(ψ::FungDemiray)
    return (:μ, :b)
end

"""
Vito [^1]

Parameters: α, β, γ

Model: ``\\alpha (\\exp\\bigg(\\beta (I_1 - 3)\\bigg) + \\gamma  (I_2 - 3)) - 1)``

[^1]: > Vito R. A note on arterial elasticity. Journal of Biomechanics. 1973 Sep 1;6(5):561-4.
"""
struct Vito end

function StrainEnergyDensityFunction(ψ::Vito, (; α, β, γ))
    W(λ⃗) = α * (exp(β * (I₁(λ⃗) - 3) + γ * (I₂(λ⃗) - 3)) - 1)
end

function StrainEnergyDensityFunction(ψ::Vito, (; α, β, γ))
    W(I⃗) = α * (exp(β * (I⃗[1] - 3) + γ * (I⃗[2] - 3)) - 1)
end

function parameters(ψ::Vito)
    return (:α, :β, :γ)
end

"""
Modified Yeoh [^1]

Parameters: C10, C20, C30, α, β

Model: ``C_{10} * (I_1 - 3) + C_{20} * (I_1 - 3)^2 + C_{30} * (I_1 - 3)^3 + \\alpha / \\beta * (1 - \\exp{-\\beta * (I_1 - 3)})``

[^1]: > He H, Zhang Q, Zhang Y, Chen J, Zhang L, Li F. A comparative study of 85 hyperelastic constitutive models for both unfilled rubber and highly filled rubber nanocomposite material. Nano Materials Science. 2021 Jul 16.
"""
struct ModifiedYeoh end

function StrainEnergyDensityFunction(ψ::ModifiedYeoh, (; C10, C20, C30, α, β))
    W(λ⃗) = C10 * (I₁(λ⃗) - 3) + C20 * (I₁(λ⃗) - 3)^2 + C30 * (I₁(λ⃗) - 3)^3 + α / β * (1 - exp(-β * (I₁(λ⃗) - 3)))
end

function StrainEnergyDensityFunction(ψ::ModifiedYeoh, (; C10, C20, C30, α, β))
    W(I⃗) = C10 * (I⃗[1] - 3) + C20 * (I⃗[1] - 3)^2 + C30 * (I⃗[1] - 3)^3 + α / β * (1 - exp(-β * (I⃗[1] - 3)))
end

function parameters(ψ::ModifiedYeoh)
    return (:C10, :C20, :C30, :α, :β)
end

"""
Chevalier-Marco [^1]

Parameters: aᵢ, bᵢ

Model: ``W = \\int\\limits_{3}^{I_1(\\vec\\lambda)} \\exp\\bigg(\\sum\\limits_{i=0}^{N}a_i(I_1-3)^i\\bigg)\\text{d}I_1+ \\int\\limits_{3}^{I_2(\\vec\\lambda)} \\sum\\limits_{i=0}^{n}\\frac{b_i}{I_2^i}\\text{d}I_2``

* NOTE: This model is not yet compatible with AD. Use Finite differences to calculate the derivatives.

[^1]: > Chevalier L, Marco Y. Tools for multiaxial validation of behavior laws chosen for modeling hyper‐elasticity of rubber‐like materials. Polymer Engineering & Science. 2002 Feb;42(2):280-98.
"""
struct ChevalierMarco end

function StrainEnergyDensityFunction(ψ::ChevalierMarco, (; aᵢ, bᵢ))
    ∂W∂I1(I₁) = exp(sum(@tullio _ := aᵢ[i] * (I₁ - 3)^(i - 1)))
    ∂W∂I2(I₂) = @tullio _ := b[i] / I₂^(i - 1)
    W(λ⃗) = quadgk(∂W∂I1, 3, I₁(λ⃗))[1] + quadgk(∂W∂I2, 3, I₂(λ⃗))[1]
end

function StrainEnergyDensityFunction(ψ::ChevalierMarco, (; aᵢ, bᵢ), I::InvariantForm)
    ∂W∂I1(I₁) = exp(sum(@tullio _ := aᵢ[i] * (I₁ - 3)^(i - 1)))
    ∂W∂I2(I₂) = @tullio _ := b[i] / I₂^(i - 1)
    W(I⃗) = quadgk(∂W∂I1, 3, I⃗[1])[1] + quadgk(∂W∂I2, 3, I⃗[2])[1]
end

function NominalStressFunction(ψ::ChevalierMarco, (; aᵢ, bᵢ))
    ∂W∂I1(λ⃗) = exp(sum(@tullio _ := aᵢ[i] * (I₁(λ⃗) - 3)^(i - 1)))
    ∂W∂I2(λ⃗) = @tullio _ := b[i] / I₂(λ⃗)^(i - 1)
    function s(λ⃗)
        sᵢ = map(λ⃗ᵢ -> 2 * (I * ∂W∂I1 - diagm(λ⃗ .^ 2)^(-2) * ∂W∂I2), λ⃗)
        sᵢ = sᵢ .- getindex.(sᵢ, 3) .* getindex.(λ⃗, 3) .* map(λ⃗ᵢ -> 1 ./ λ⃗ᵢ, λ⃗)
        return sᵢ
    end
end

function TrueStressFunction(ψ::ChevalierMarco, (; aᵢ, bᵢ))
    ∂W∂I1(λ⃗) = exp(sum(@tullio _ := aᵢ[i] * (I₁(λ⃗) - 3)^(i - 1)))
    ∂W∂I2(λ⃗) = @tullio _ := b[i] / I₂(λ⃗)^(i - 1)
    s = NominalStressFunction(ψ, (aᵢ=aᵢ, bᵢ=bᵢ))
    function σ(λ⃗)
        σᵢ = map(λ⃗ᵢ -> λ⃗ᵢ .* s(λ⃗ᵢ), λ⃗)
        return σᵢ
    end
end

function parameters(ψ::ChevalierMarco)
    return (:aᵢ, :bᵢ)
end

"""
Gornet - Desmorat [^1]

Parameters: h₁, h₂, h₃

Model: ``W = h_1\\int\\exp{h_3(I_1-3)^2}\\text{d}I_1+3h_2\\int\\frac{1}{\\sqrt{I_2}}\\text{d}I_2 = \\frac{h_1 \\sqrt{\\pi} \\text{erfi}(\\sqrt{h_3}(I_1-3)^2)}{2\\sqrt{h_3}}+6h_2\\sqrt{I_2}``

* Note: the differential form was original form and the closed form SEF was determine via symbolic integration in Mathematica.

[^1]: > Gornet L, Marckmann G, Desmorat R, Charrier P. A new isotropic hyperelastic strain energy function in terms of invariants and its derivation into a pseudo-elastic model for Mullins effect: application to finite element analysis. Constitutive Models for Rubbers VII. 2012:265-71.
"""
struct GornetDesmorat end

function StrainEnergyDensityFunction(ψ::GornetDesmorat, (; h₁, h₂, h₃))
    W(λ⃗) = h₁ * √π * erfi(√h₃ * (I₁(λ⃗) - 3)^2) / 2 / √h₃ + 6 * h₂ * √(I₂(λ⃗))
end

function StrainEnergyDensityFunction(ψ::GornetDesmorat, (; h₁, h₂, h₃), I::InvariantForm)
    W(I⃗) = h₁ * √π * erfi(√h₃ * (I⃗[1] - 3)^2) / 2 / √h₃ + 6 * h₂ * √(I⃗[2])
end

function parameters(ψ::GornetDesmorat)
    return (:h₁, :h₂, :h₃)
end

"""
Mansouri-Darijani [^1]

Parameters: A1, m1, B1, n1

Model: ``A_1\\exp{m_1(I_1-3)-1}+B_1\\exp{n_1(I_2-3)-1}``

[^1]: > Mansouri MR, Darijani H. Constitutive modeling of isotropic hyperelastic materials in an exponential framework using a self-contained approach. International Journal of Solids and Structures. 2014 Dec 1;51(25-26):4316-26.
"""
struct MansouriDarijani end

function StrainEnergyDensityFunction(ψ::MansouriDarijani, (; A1, m1, B1, n1))
    W(λ⃗) = A1 * (exp(m1 * (I₁(λ⃗) - 3)) - 1) + B1 * (exp(n1 * (I₂(λ⃗) - 3)) - 1)
end

function StrainEnergyDensityFunction(ψ::MansouriDarijani, (; A1, m1, B1, n1), I::InvariantForm)
    W(I⃗) = A1 * (exp(m1 * (I⃗[1] - 3)) - 1) + B1 * (exp(n1 * (I⃗[2] - 3)) - 1)
end

function parameters(ψ::MansouriDarijani)
    return (:A1, :m1, :B1, :n1)
end

"""
Gent Thomas [^1]

Paramters: C1, C2

Model: ``C_1(I_1-3)+C_2\\log(\\frac{I_2}{3})``

[^1]: > Gent AN, Thomas AG. Forms for the stored (strain) energy function for vulcanized rubber. Journal of Polymer Science. 1958 Apr;28(118):625-8.
"""
struct GentThomas end

function StrainEnergyDensityFunction(ψ::GentThomas, (; C1, C2))
    W(λ⃗) = C1 * (I₁(λ⃗) - 3) + C2 * log(I₂(λ⃗) / 3)
end

function StrainEnergyDensityFunction(ψ::GentThomas, (; C1, C2), I::InvariantForm)
    W(I⃗) = C1 * (I⃗[1] - 3) + C2 * log(I⃗[2] / 3)
end

function parameters(ψ::GentThomas)
    return (:C1, :C2)
end

"""
Alexander [^1]

Parameters: C₁, C₂, C₃, k, γ

Model: ``\\frac{C_1 \\sqrt{\\pi}\\text{erfi}\\big(\\sqrt{k}(I_1-3)\\big)}{2\\sqrt{k}}+C_2\\log{\\frac{I_2-3+\\gamma}{\\gamma}}+C_3(I_2-3)``

[^1]: > Alexander H. A constitutive relation for rubber-like materials. International Journal of Engineering Science. 1968 Sep 1;6(9):549-63.
"""
struct Alexander end

function StrainEnergyDensityFunction(ψ::Alexander, (; C₁, C₂, C₃, k, γ))
    W(λ⃗) = C₁ * √π * erfi(√k * (I₁(λ⃗) - 3)) / 2 / √k + C₂ * log((I₂(λ⃗) - 3 + γ) / γ) + C₃ * (I₂(λ⃗) - 3)
end

function StrainEnergyDensityFunction(ψ::Alexander, (; C₁, C₂, C₃, k, γ), I::InvariantForm)
    W(I⃗) = C₁ * √π * erfi(√k * (I⃗[1] - 3)) / 2 / √k + C₂ * log((I⃗[2] - 3 + γ) / γ) + C₃ * (I⃗[2] - 3)
end

function parameters(ψ::Alexander)
    return (:C₁, :C₂, :C₃, :k, :γ)
end

"""
Lambert-Diani Rey [^1]

Parameters: aᵢ, bᵢ

Model: ``\\int\\limits_{3}^{I_1}\\exp\\bigg(\\sum\\limits_{i=0}^{n}a_i(I_1-3)^i\\bigg)\\text{d}I_1+\\int\\limits_{3}^{I_2}\\sum\\limits_{j=0}^{m}b_i\\log(I_2)^i\\text{d}I_2``

[^1]: > Lambert-Diani J, Rey C. New phenomenological behavior laws for rubbers and thermoplastic elastomers. European Journal of Mechanics-A/Solids. 1999 Nov 1;18(6):1027-43.
"""
struct LambertDianiRey end

function StrainEnergyDensityFunction(ψ::LambertDianiRey, (; aᵢ, bᵢ))
    function W(λ⃗)
        ∂W∂I₁(I₁) = exp(sum(aᵢ .* (I₁ .- 3) .^ i))
        ∂W∂I₂(I₂) = exp(sum(b₁ .* log(I₂) .^ i))
        W(λ⃗) = quadgk(∂W∂I₁, 3, I₁(λ⃗))[1] + quadgk(∂W∂I₂, 3, I₂(λ⃗))[1]
    end
    return W
end

function StrainEnergyDensityFunction(ψ::LambertDianiRey, (; aᵢ, bᵢ), I::InvariantForm)
    function W(I⃗)
        ∂W∂I₁(I₁) = exp(sum(aᵢ .* (I₁ .- 3) .^ i))
        ∂W∂I₂(I₂) = exp(sum(b₁ .* log(I₂) .^ i))
        W(I⃗) = quadgk(∂W∂I₁, 3, I⃗[1])[1] + quadgk(∂W∂I₂, 3, I⃗[2])[1]
    end
    return W
end


function NominalStressFunction(ψ::LambertDianiRey, (; aᵢ, bᵢ))
    function s(λ⃗)
        ∂W∂I₁ = exp(sum(aᵢ .* (I₁(λ⃗) .- 3) .^ i))
        ∂W∂I₂ = exp(sum(b₁ .* log(I₂(λ⃗)) .^ i))
        sᵢ = map(λ⃗ᵢ -> 2 * (I * ∂W∂I₁ - diagm(λ⃗ .^ 2)^(-2) * ∂W∂I₂), λ⃗)
        sᵢ = sᵢ .- getindex.(sᵢ, 3) .* getindex.(λ⃗, 3) .* map(λ⃗ᵢ -> 1 ./ λ⃗ᵢ, λ⃗)
        return sᵢ
    end
    return s
end

function TrueStressFunction(ψ::LambertDianiRey, (; aᵢ, bᵢ))
    s = NominalStressFunction(ψ, (aᵢ=aᵢ, bᵢ=bᵢ))
    function σ(λ⃗)
        σᵢ = map(λ⃗ᵢ -> λ⃗ᵢ .* s(λ⃗ᵢ), λ⃗)
        return σᵢ
    end
    return σᵢ
end

function parameters(ψ::LambertDianiRey)
    return (:aᵢ, :bᵢ)
end

"""
Hoss Marczak I [^1]

Parameters: α, β, μ, b, n

Model: ``\\frac{\\alpha}{\\beta}(1-\\exp{-\\beta(I_1-3)})+\\frac{\\mu}{2b}\\bigg((1+\\frac{b}{n}(I_1-3))^n -1\\bigg)``

* Note: The authors suggested this model for low strains.

[^1]: > Hoss L, Marczak RJ. A new constitutive model for rubber-like materials. Mecánica Computacional. 2010;29(28):2759-73.
"""
struct HossMarczakI end

function StrainEnergyDensityFunction(ψ::HossMarczakI, (; α, β, μ, b, n))
    W(λ⃗) = α / β * (1 - exp(-β * (I₁(λ⃗) - 3))) + μ / (2b) * ((1 + b / n * (I₁(λ⃗) - 3))^n - 1)
end

function StrainEnergyDensityFunction(ψ::HossMarczakI, (; α, β, μ, b, n), I::InvariantForm)
    W(I⃗) = α / β * (1 - exp(-β * (I⃗[1] - 3))) + μ / (2b) * ((1 + b / n * (I⃗[1] - 3))^n - 1)
end

function parameters(ψ::HossMarczakI)
    return (:α, :β, :μ, :b, :n)
end

"""
Hoss Marczak II [^1]

Parameters: α, β, μ, b, n, C2

Model: ``\\frac{\\alpha}{\\beta}(1-\\exp{-\\beta(I_1-3)})+\\frac{\\mu}{2b}\\bigg((1+\\frac{b}{n}(I_1-3))^n -1\\bigg)+C_2\\log(\\frac{I_2}{3})``

* Note: The authors suggests this model for high strains.

[^1]: > Hoss L, Marczak RJ. A new constitutive model for rubber-like materials. Mecánica Computacional. 2010;29(28):2759-73.
"""
struct HossMarczakII end

function StrainEnergyDensityFunction(ψ::HossMarczakII, (; α, β, μ, b, n, C2))
    W(λ⃗) = α / β * (1 - exp(-β * (I₁(λ⃗) - 3))) + μ / (2b) * ((1 + b / n * (I₁(λ⃗) - 3))^n - 1) + C2 * log(I₂(λ⃗) / 3)
end

function StrainEnergyDensityFunction(ψ::HossMarczakII, (; α, β, μ, b, n, C2), I::InvariantForm)
    W(I⃗) = α / β * (1 - exp(-β * (I⃗[1] - 3))) + μ / (2b) * ((1 + b / n * (I⃗[1] - 3))^n - 1) + C2 * log(I⃗[2] / 3)
end

function parameters(ψ::HossMarczakII)
    return (:α, :β, :μ, :b, :n, :C2)
end

"""
Exp-Ln [^1]

Parameters: A, a, b

Model: ``A\\bigg[\\frac{1}{a}\\exp{(a(I_1-3))}+b(I_1-2)(1-\\log{I_1-2})-\\frac{1}{a}-b\\bigg]``

[^1]: > Khajehsaeid H, Arghavani J, Naghdabadi R. A hyperelastic constitutive model for rubber-like materials. European Journal of Mechanics-A/Solids. 2013 Mar 1;38:144-51.
"""
struct ExpLn end

function StrainEnergyDensityFunction(ψ::ExpLn, (; A, a, b))
    W(λ⃗) = A * (1 / a * exp(a * (I₁(λ⃗) - 3)) + b * (I₁(λ⃗) - 2) * (1 - log(I₁(λ⃗) - 2)) - 1 / a - b)
end

function StrainEnergyDensityFunction(ψ::ExpLn, (; A, a, b), I::InvariantForm)
    W(I⃗) = A * (1 / a * exp(a * (I⃗[1] - 3)) + b * (I⃗[1] - 2) * (1 - log(I⃗[1] - 2)) - 1 / a - b)
end

function parameters(ψ::ExpLn)
    return (:A, :a, :b)
end

"""
Van der Waals [^1][^2][^3]

Parameters: μ, λm, β, α

Model:

``W(\\vec{\\lambda}) = -\\mu\\{(\\lambda_m^2-3)\\log(1-\\Theta)+\\Theta\\}-\\frac{2\\alpha}{3}\\bigg(\\frac{I-3}{2}\\bigg)^{3/2}``

``\\Theta = \\frac{\\beta I_1 + (1-\\beta)I_2-3}{\\lambda_m^2-3)}``

[^1]: > Kilian HG, Enderle HF, Unseld K. The use of the van der Waals model to elucidate universal aspects of structure-property relationships in simply extended dry and swollen rubbers. Colloid and Polymer Science. 1986 Oct;264(10):866-76.
[^2]: > Ambacher H, Enderle HF, Kilian HG, Sauter A. Relaxation in permanent networks. InRelaxation in Polymers 1989 (pp. 209-220). Steinkopff.
[^3]: > Kilian HG. A molecular interpretation of the parameters of the van der Waals equation of state for real networks. Polymer Bulletin. 1980 Sep;3(3):151-8.
"""
struct VanDerWaals end

function StrainEnergyDensityFunction(ψ::VanDerWaals, (; μ, λm, β, α))
    function W(λ⃗)
        I = β * I₁(λ⃗) + (1 - β) * I₂(λ⃗)
        θ = (I - 3) / (λm^2 - 3)
        μ * (-(λm^2 - 3) * log(1 - θ) + θ) - 2 / 3 * α * ((I - 3) / 2)^(3 / 2)
    end
end

function StrainEnergyDensityFunction(ψ::VanDerWaals, (; μ, λm, β, α), I::InvariantForm)
    function W(I⃗)
        I = β * I⃗[1] + (1 - β) * I⃗[2]
        θ = (I - 3) / (λm^2 - 3)
        μ * (-(λm^2 - 3) * log(1 - θ) + θ) - 2 / 3 * α * ((I - 3) / 2)^(3 / 2)
    end
end


function parameters(ψ::VanDerWaals)
    return (:μ, :λm, :β, :α)
end

"""
Gent [^1]

Parameters: μ, Jₘ

Model: ``-\\frac{\\mu J_m}{2}\\log{\\bigg(1-\\frac{I_1-3}{J_m}\\bigg)}``

[^1]: > Gent AN. A new constitutive relation for rubber. Rubber chemistry and technology. 1996 Mar;69(1):59-61.
"""
struct Gent end

function StrainEnergyDensityFunction(ψ::Gent, (; μ, Jₘ))
    W(λ⃗) = -(μ * Jₘ) / 2 * log(1 - (I₁(λ⃗) - 3) / Jₘ)
end

function StrainEnergyDensityFunction(ψ::Gent, (; μ, Jₘ), I::InvariantForm)
    W(I⃗) = -(μ * Jₘ) / 2 * log(1 - (I⃗[1] - 3) / Jₘ)
end

function parameters(ψ::Gent)
    return (:μ, :Jₘ)
end

"""
Takamizawa-Hayashi [^1]
From: A description of arterial wall mechanics using limiting chain extensibility constitutitive models by Horgan and Saccomandi

Parameters: c, Jₘ

Model: ``-c\\log{1-\\big(\\frac{I_1-3}{J_m}\\big)^2}``

[^1]: > Takamizawa K, Hayashi K. Strain energy density function and uniform strain hypothesis for arterial mechanics. Journal of biomechanics. 1987 Jan 1;20(1):7-17.
"""
struct TakamizawaHayashi end

function StrainEnergyDensityFunction(ψ::TakamizawaHayashi, (; c, Jₘ))
    W(λ⃗) = -c * log(1 - ((I₁(λ⃗) - 3) / Jₘ)^2)
end

function StrainEnergyDensityFunction(ψ::TakamizawaHayashi, (; c, Jₘ), I::InvariantForm)
    W(I⃗) = -c * log(1 - ((I⃗[1] - 3) / Jₘ)^2)
end

function parameters(ψ::TakamizawaHayashi)
    return (:c, :Jₘ)
end

"""
Yeoh-Fleming [^1]

Parameters: A, B, C10, Im

Model: ``\\frac{A}{B}(1-\\exp{-B(I_1-3)}) - C_{10}(I_m-3)\\log{1-\\frac{I_1-3}{I_m-3}}``

[^1]: >  Yeoh OH, Fleming PD. A new attempt to reconcile the statistical and phenomenological theories of rubber elasticity. Journal of Polymer Science Part B: Polymer Physics. 1997 Sep 15;35(12):1919-31.
"""
struct YeohFleming end

function StrainEnergyDensityFunction(ψ::YeohFleming, (; A, B, C10, Im))
    W(λ⃗) = A / B * (1 - exp(-B * (I₁(λ⃗) - 3))) - C10 * (Im - 3) * log(1 - ((I₁(λ⃗) - 3) / (Im - 3)))
end

function StrainEnergyDensityFunction(ψ::YeohFleming, (; A, B, C10, Im), I::InvariantForm)
    W(I⃗) = A / B * (1 - exp(-B * (I⃗[1] - 3))) - C10 * (Im - 3) * log(1 - ((I⃗[1] - 3) / (Im - 3)))
end

function parameters(ψ::YeohFleming)
    return (:A, :B, :C10, :Im)
end

"""
Pucci-Saccomandi [^1]

Parameters: K, μ, Jₘ

Model ``K\\log{\\frac{I_2}{3}}-\\frac{\\mu J_m}{2}\\log{1-\\frac{I_1-3}{J-m}}``

[^1]: > Pucci E, Saccomandi G. A note on the Gent model for rubber-like materials. Rubber chemistry and technology. 2002 Nov;75(5):839-52.
"""
struct PucciSaccomandi end

function StrainEnergyDensityFunction(ψ::PucciSaccomandi, (; K, μ, Jₘ))
    W(λ⃗) = K * log(I₂(λ⃗) / 3) - μ * Jₘ / 2 * log(1 - (I₁(λ⃗) - 3) / Jₘ)
end

function StrainEnergyDensityFunction(ψ::PucciSaccomandi, (; K, μ, Jₘ), I::InvariantForm)
    W(I⃗) = K * log(I⃗[2] / 3) - μ * Jₘ / 2 * log(1 - (I⃗[1] - 3) / Jₘ)
end

function parameters(ψ::PucciSaccomandi)
    return (:K, :μ, :Jₘ)
end

# Originally from CONSTITUTIVE MODELS FOR ATACTIC ELASTOMERS
"""
Horgan Saccomandi Model

Parameters: μ, J

Model: ``-\\frac{\\mu J}{2}\\log\\bigg(\\frac{J^3-J^2I_1+JI_2-1}{(J-1)^3}\\bigg)``
"""
struct HorganSaccomandi end

function StrainEnergyDensityFunction(ψ::HorganSaccomandi, (; μ, J))
    W(λ⃗) = -μ * J / 2 * log((J^3 - J^2 * I₁(λ⃗) + J * I₂(λ⃗) - 1) / (J - 1)^3)
end

function StrainEnergyDensityFunction(ψ::HorganSaccomandi, (; μ, J), I::InvariantForm)
    W(I⃗) = -μ * J / 2 * log((J^3 - J^2 * I⃗[1] + J * I⃗[2] - 1) / (J - 1)^3)
end

function parameters(ψ::HorganSaccomandi)
    return (:μ, :J)
end

"""
Beatty Model

Parameters: G₀, Iₘ

Model: ``-\\frac{G_0 I_m(I_m-3)}{2(2I_m-3)}\\log\\bigg(\\frac{1-\\frac{I_1-3}{I_m-3}}{1+\\frac{I_1-3}{I_m}} \\bigg)``
"""
struct Beatty end

function StrainEnergyDensityFunction(ψ::Beatty, (; G₀, Iₘ))
    W(λ⃗) = -G₀ * Iₘ * (Iₘ - 3) / 2 / (2Iₘ - 3) * log((1 - (I₁(λ⃗) - 3) / (Iₘ - 3)) / (1 + (I₁(λ⃗) - 3) / (Iₘ)))
end

function StrainEnergyDensityFunction(ψ::Beatty, (; G₀, Iₘ), I::InvariantForm)
    W(I⃗) = -G₀ * Iₘ * (Iₘ - 3) / 2 / (2Iₘ - 3) * log((1 - (I⃗[1] - 3) / (Iₘ - 3)) / (1 + (I⃗[1] - 3) / (Iₘ)))
end

function parameters(ψ::Beatty)
    return (:G₀, :Iₘ)
end

"""
Horgan Murphy Model

Parameters: μ, Jₘ, c

Model: ``-\\frac{2\\mu J_m}{c^2}\\log\\bigg(1-\\frac{\\lambda_1^c+\\lambda_2^c+\\lambda_3^c-3}{J_m})``
"""
struct HorganMurphy end

function StrainEnergyDensityFunction(ψ::HorganMurphy, (; μ, Jₘ, c))
    W(λ⃗) = -2 * μ * Jₘ / c^2 * log(1 - (sum(λ⃗ .^ c) - 3) / Jₘ)
end

function parameters(ψ::HorganMurphy)
    return (:μ, :Jₘ, :c)
end

########################
########### TABLE 4
########################
"""
Valanis-Landel

Parameters: μ

Model: ``2\\mu\\sum\\limits_{1}^{3}(\\lambda_i(\\log\\lambda_i -1))``
"""
struct ValanisLandel end

function StrainEnergyDensityFunction(ψ::ValanisLandel, (; μ))
    W(λ⃗) = 2 * μ * sum(λ⃗ * (log.(λ⃗) - 1))
end

function parameters(ψ::ValanisLandel)
    return (:μ)
end

"""
Peng - Landel

Parameters: E

Model: ``E\\sum\\limits_{i=1}^{3}\\bigg[\\lambda_i - 1 - \\log(\\lambda_i) - \\frac{1}{6}\\log(\\lambda_i)^2 + \\frac{1}{18}\\log(\\lambda_i)^3-\\frac{1}{216}\\log(\\lambda_i)^4\\bigg]``
"""
struct PengLandel end

function StrainEnergyDensityFunction(ψ::PengLandel, (; E))
    W(λ⃗) = sum(@. λ⃗ - 1 - log(λ⃗) - 1 / 6 * log(λ⃗)^2 + 1 / 18 * log(λ⃗)^3 - 1 / 216 * log(λ⃗)^4) * E
end

function parameters(ψ::PengLandel)
    return (:E)
end


"""
Ogden

Parameters: μ⃗, α⃗

Model: ``\\sum\\limits_{i=1}^{N}\\frac{\\mu_i}{\\alpha_i}(\\lambda_1^{\\alpha_i}+\\lambda_2^{\\alpha_i}+\\lambda_3^{\\alpha_i}-3)``
"""
struct Ogden end

function StrainEnergyDensityFunction(ψ::Ogden, (; μ, α))
    W(λ⃗) = @tullio _ := μ[i] / α[i] * (sum(λ⃗ .^ α[i]) - 3)
end

function parameters(ψ::Ogden)
    return (:μ, :α)
end

"""
Attard

Parameters: A⃗, B⃗

Model: ``\\sum\\limits_{i=1}^N\\frac{A_i}{2i}(\\lambda_1^{2i}+\\lambda_2^{2i}+\\lambda_3^{2i}-3) + \\frac{B_i}{2i}(\\lambda_1^{-2i}+\\lambda_2^{-2i}+\\lambda_3^{-2i}-3)``
"""
struct Attard end

function StrainEnergyDensityFunction(ψ::Attard, (; A, B))
    @assert length(A)==length(B) "Length of A and B are not equal"
    W(λ⃗) = @tullio _ := A[i] / 2 / i * (sum(λ⃗ .^ (2i)) - 3) + B[i] / 2 / i * (sum(λ⃗ .^ (-2i)) - 3)
end

function parameters(ψ::Attard)
    return (:A, :B)
end

"""
Shariff

Parameters: E, α₁, α₂, α₃, α₄, α₅

Model:
``E\\sum\\limits_{i=1}^3\\sum\\limits_{j=1}^{N}\\alpha_j \\Phi_j(\\lambda_i)``
"""
struct Shariff end

function StrainEnergyDensityFunction(ψ::Shariff, (; E, α))
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

function parameters(ψ::Shariff)
    return (:E, :α)
end


# Article requested
"""
Arman - Narooei

Parameters: A⃗, B⃗, m⃗, n⃗, α⃗, β⃗

Model: ``\\sum\\limits_{i=1}^{N} A_i\\big[\\exp{m_i(\\lambda_1^{\\alpha_i}+\\lambda_2^{\\alpha_i}+\\lambda_3^{\\alpha_i}-3)}-1] + B_i\\big[\\exp{n_i(\\lambda_1^{-\\beta_i}+\\lambda_2^{-\\beta_i}+\\lambda_3^{-\\beta_i}-3)}-1]``
"""
struct ArmanNarooei end

function StrainEnergyDensityFunction(ψ::ArmanNarooei, (; A, B, m, n, α, β))
    @assert length(A)==length(B)==length(m)==length(n)==length(α)==length(β) "Length of A, B, m, n, α and β are not equal"
    W(λ⃗) = @tullio _ := A[i] * (exp(m[i] * (sum(λ⃗ .^ α[i]) - 3)) - 1) + B[i] * (exp(n[i] * (sum(λ⃗ .^ (-β[i])) - 3)) - 1)
end

function parameters(ψ::ArmanNarooei)
    return (:A, :B, :m, :n, :α, :β)
end

################################
###################      Table 5
################################
"""
Continuum Hybrid

Parameters: K₁, K₂, α, μ

Model: ``K_1(I_1-3)+K_2\\log\\frac{I_2}{3}+\\frac{\\mu}{\\alpha}(\\lambda_1^\\alpha+\\lambda_2^\\alpha+\\lambda^\\alpha-3)``
"""
struct ContinuumHybrid end

function StrainEnergyDensityFunction(ψ::ContinuumHybrid, (; K₁, K₂, α, μ))
    W(λ⃗) = K₁ * (I₁(λ⃗) - 3) + K₂ * log(I₂(λ⃗) / 3) + μ / α * (sum(λ⃗ .^ α) - 3)
end

function parameters(ψ::ContinuumHybrid)
    return (:K₁, :K₂, :α, :μ)
end

"""
Bechir-4 Term

Parameters: C11, C12, C21, C22

Model: ``C_1^1(I_1-3)+\\sum\\limits_{n=1}^{2}\\sum\\limits_{r=1}^{2}C_n^{r}(\\lambda_1^{2n}+\\lambda_2^{2n}+\\lambda_3^{2n}-3)^r``
"""
struct Bechir4Term end

function StrainEnergyDensityFunction(ψ::Bechir4Term, (; C11, C12, C21, C22))
    C = [C11 C12; C21 C22]
    W(λ⃗) = C[1, 1] * (I₁(λ⃗) - 3) + sum(n -> sum(r -> C[n, r] * (sum(λ⃗ .^ (2n))), 1:2), 1:2)
end

function parameters(ψ::Bechir4Term)
    return (:C11, :C12, :C21, :C22)
end

"""
Constrained Junction

Parameters: Gc, νkT, κ

Model: ``G_c (I_1-3)+ \\frac{\\nu k T}{2}(\\sum\\limits_{i=1}^{3}\\kappa\\frac{\\lambda_i-1}{\\lambda_i^2+\\kappa}+\\log{\\frac{\\lambda_i^2+\\kappa}{1+\\kappa}}-\\log{\\lambda_i^2})``
"""
struct ConstrainedJunction end

function StrainEnergyDensityFunction(ψ::ConstrainedJunction, (; Gc, νkT, κ))
    W(λ⃗) = Gc * (I₁(λ⃗) - 3) + μkT / 2 * sum(i -> κ * (λ⃗[i] - 1) / (λ⃗[i]^2 + κ) + log((λ⃗[i]^2 + κ) / (1 + κ)) - log(λ⃗[i]^2), 1:3)
end

function parameters(ψ::ConstrainedJunction)
    return (:Gc, :νkT, :κ)
end

"""
Edward-Vilgis

Parameters: Ns, Nc, α, η

Model: ``\\frac{1}{2}N_C\\Bigg[\\frac{(1-\\alpha^2)I_1}{1-\\alpha^2I_1}+\\log(1-\\alpha^2I_1)\\Bigg]+\\frac{1}{2}N_S\\Bigg[\\sum_{i=1}^{3}\\Big\\{\\frac{(1+\\eta)(1-\\alpha^2)\\lambda_i^2}{( 1+\\eta\\lambda_i^2)(1-\\alpha^2I_1)}+\\log(1+\\eta\\lambda_i^2)\\Big\\}+\\log(1-\\alpha^2I_1)\\Bigg]``
"""
struct EdwardVilgis end

function StrainEnergyDensityFunction(ψ::EdwardVilgis, (; Ns, Nc, α, η))
    W(λ⃗) = 0.5 * Nc * ((1 - α^2) * I₁(λ⃗) / (1 - α^2 * I₁(λ⃗)) + log(1 - α^2 * I₁(λ⃗))) + 0.5 * Ns * ((1 + η) * (1 - α^2) * λ⃗[1] / (1 + η * λ⃗[1]^2) / (1 - α^2 * I₁(λ⃗)) + log(1 + η * λ⃗[1]^2) + (1 + η) * (1 - α^2) * λ⃗[2] / (1 + η * λ⃗[2]^2) / (1 - α^2 * I₁(λ⃗)) + log(1 + η * λ⃗[2]^2) + (1 + η) * (1 - α^2) * λ⃗[3] / (1 + η * λ⃗[3]^2) / (1 - α^2 * I₁(λ⃗)) + log(1 + η * λ⃗[3]^2) + log(1 - α^2 * I₁(λ⃗)))
end

function parameters(ψ::EdwardVilgis)
    return (:Ns, :Nc, :α, :η)
end

"""
MCC (modified constrained chain)

Parameters:

Model:

``\\frac{1}{2}\\zeta k T \\sum\\limits_{i=1}^{3}(\\lambda_i^2-1)+\\frac{1}{2}\\mu k T\\sum\\limits_{i=1}^{3}[B_i+D_i-\\log{(1+B_i)}-\\log{(1+D_i)}]``

``B_i = \\frac{\\kappa^2(\\lambda_i^2-1)}{(\\lambda_i^2+\\kappa)^2}``

``D_i = \\frac{\\lambda_i^2 B_i}{\\kappa}``
"""
struct MCC end

function StrainEnergyDensityFunction(ψ::MCC, (; ζkT, μkT, κ))
    W(λ⃗) = 1 / 2 * ζkT * sum(i -> λ⃗[i]^2 - 1, 1:3) + 1 / 2 * μkT * sum(i -> κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2) + (λ⃗[i]^2 * (κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2)) / κ) - log(1 + (κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2))) - log(1 + (λ⃗[i]^2 * (κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2)) / κ)))
end

function parameters(ψ::MCC)
    return (:ζkT, :μkT, :κ)
end

"""
Tube

Parameters: Gc, Ge, β

Model: ``\\sum\\limits_{i=1}^{3}\\frac{G_c}{2}(\\lambda_i^2-1)+\\frac{2Ge}{\\beta^2}(\\lambda_i^{-\\beta}-1)``
"""
struct Tube end

function StrainEnergyDensityFunction(ψ::Tube, (; Gc, Ge, β))
    W(λ⃗) = @tullio _ := Gc / 2 * (λ⃗[i]^2 - 1) + 2Ge / β^2 * (λ⃗[i]^(-β) - 1)
end

function parameters(ψ::Tube)
    return (:Gc, :Ge, :β)
end

"""
Nonaffine - Tube

Parameters: Gc, Ge

Model: ``G_c \\sum\\limits_{i=1}^{3}\\frac{\\lambda_i^2}{2}+G_e\\sum\\limits_{i=1}^{3}\\lambda_i+\\frac{1}{\\lambda_i}``
"""
struct NonaffineTube end

function StrainEnergyDensityFunction(ψ::NonaffineTube, (; Gc, Ge))
    W(λ⃗) = Gc * sum(λ⃗ .^ 2 ./ 2) + Ge * sum(λ⃗ .+ 1 ./ λ⃗)
end

function parameters(ψ::NonaffineTube)
    return (:Gc, :Ge)
end

"""
Three Chain Model

Parameters: μ, N

Model: `` \\frac{\\mu\\sqrt{N}}{3}\\sum\\limits_{i=1}^{3}\\bigg(\\lambda_i\\beta_i+\\sqrt{N}\\log\\bigg(\\frac{\\beta_i}{\\sinh \\beta_i}\\bigg)\\bigg)``

"""
struct ThreeChainModel end

function StrainEnergyDensityFunction(ψ::ThreeChainModel, (; μ, N))
    ℒinv(x) = x * (3 - 1.0651 * x^2 - 0.962245 * x^4 + 1.47353 * x^6 - 0.48953 * x^8) / (1 - x) / (1 + 1.01524 * x)
    W(λ⃗) = μ * sqrt(N) / 3 * sum(λ⃗ .* ℒinv.(λ⃗ ./ sqrt(N)) .+ sqrt(N) .* log.((ℒinv.(λ⃗ ./ sqrt(N))) ./ (sinh.(ℒinv.(λ⃗ ./ sqrt(N))))))
end

function parameters(ψ::ThreeChainModel)
    return (:μ, :N)
end

"""
Arruda Boyce

Parameters: μ, N

Model: ``\\mu\\bigg(\\frac{1}{2}(I_1-3)+\\frac{I_1^2-9}{20N}+\\frac{11(I_1^3-27)}{1050N^2}+\\frac{19(I_1^4-81)}{7000N^3}+\\frac{519(I_1^5-243)}{673750N^4}\\bigg)``
"""
struct ArrudaBoyce end

function StrainEnergyDensityFunction(ψ::ArrudaBoyce, (; μ, N))
    # ℒinv(x) = x * (3 - 1.0651 * x^2 - 0.962245 * x^4 + 1.47353 * x^6 - 0.48953 * x^8) / (1 - x) / (1 + 1.01524 * x)
    ℒinv(x) = x * (3 - x^2) / (1 - x^2)
    function W(λ⃗)
        rchain_Nl = √(I₁(λ⃗) / 3 / N)
        β = ℒinv(rchain_Nl)
        μ * N * (rchain_Nl * β + log(β / sinh(β)))
    end
end

function StrainEnergyDensityFunction(ψ::ArrudaBoyce, (; μ, N), I::InvariantForm)
    # ℒinv(x) = x * (3 - 1.0651 * x^2 - 0.962245 * x^4 + 1.47353 * x^6 - 0.48953 * x^8) / (1 - x) / (1 + 1.01524 * x)
    ℒinv(x) = x * (3 - x^2) / (1 - x^2)
    function W(I⃗)
        rchain_Nl = √(I⃗[1] / 3 / N)
        β = ℒinv(rchain_Nl)
        μ * N * (rchain_Nl * β + log(β / sinh(β)))
    end
end

function parameters(ψ::ArrudaBoyce)
    return (:μ, :N)
end

"""
Modified Flory Erman

Parameters: μ, N, κ

Model: ``W_{\\text{Arruda-Boyce}}+\\sum\\limits_{i=1}^{3}\\frac{\\mu}{2}[B_i+D_i]
"""
struct ModifiedFloryErman end

function StrainEnergyDensityFunction(ψ::ModifiedFloryErman, (; μ, N, κ))
    ArrudaBoyce((μ=μ, N=N))
    WAB = StrainEnergyDensityFunction(ArrudaBoyce(), (μ=μ, N=N))
    function W(λ⃗)
        B = map(i -> κ^2 * (λ⃗[i]^2 - 1) / (λ⃗[i]^2 + κ)^2, 1:3)
        D = map(i -> λ⃗[i]^2 * B[i] / κ, 1:3)
        WAB(λ⃗) + map(i -> B[i] + D[i] - log(B[i] + 1) - log(D[i] + 1), 1:3)
    end
end

function parameters(ψ::ModifiedFloryErman)
    return (:μ, :N, :κ)
end

"""
Extended Tube Model [^1]

Parameters: Gc, Ge, δ, β

Model: ``\\frac{G_c}{2}\\bigg[\\frac{(1-\\delta^2)(I_1-3)}{1-\\delta^2(I_1-3)}+\\log{(1-\\delta^2(I_1-3))}\\bigg]+\\frac{2G_e}{\\beta^2}\\sum\\limits_{i=1}^{3}(\\lambda_i^{-\\beta}-1)``

[^1]: > Kaliske M, Heinrich G. An extended tube-model for rubber elasticity: statistical-mechanical theory and finite element implementation. Rubber Chemistry and Technology. 1999 Sep;72(4):602-32.
"""
struct ExtendedTubeModel end

function StrainEnergyDensityFunction(ψ::ExtendedTubeModel, (; Gc, Ge, δ, β))
    W(λ⃗) = Gc / 2 * ((1 - δ^2) * (I₁(λ⃗) - 3) / (1 - δ^2 * (I₁ - 3)) + log(1 - δ^2 * (I₁(λ⃗) - 3))) + 2 * Ge / β^2 * sum(λ⃗ .^ (-β) .- 1)
end

function parameters(ψ::ExtendedTubeModel)
    return (:Gc, :Ge, :δ, :β)
end

"""
ABGI

Parameters: μ, N, Ge, n

Model: ``W_{Arruda-Boyce} + G_e\\frac{\\lambda_1^n+\\lambda_2^2+\\lambda_3^2-3}{n}``
"""
struct ABGI end

function StrainEnergyDensityFunction(ψ::ABGI, (; μ, N, Ge, n))
    WAB = StrainEnergyDensityFunction(ArrudaBoyce(), (μ=μ, N=N))
    W(λ⃗) = WAB(λ⃗) + Ge * (sum(λ⃗ .^ n) - 3) / n
end

function parameters(ψ::ABGI)
    return (:μ, :N, :Ge, :n)
end

"""
Non-Affine Micro-Sphere [^1]

Parameters: μ, N, p, U, q

Model: See Paper

---
[^1]: > Miehe C, Göktepe S, Lulei F. A micro-macro approach to rubber-like materials—part I: the non-affine micro-sphere model of rubber elasticity. Journal of the Mechanics and Physics of Solids. 2004 Nov 1;52(11):2617-60.
"""
struct NonaffineMicroSphere end

function StrainEnergyDensityFunction(ψ::NonaffineMicroSphere, (; μ, N, p, U, q))
    a = √(2) / 2
    b = 0.836095596749
    c = 0.387907304067
    r⃗ = [
        [0, 0, 1],
        [0, 1, 0],
        [1, 0, 0],
        [0, a, a],
        [0, -a, a],
        [a, 0, a],
        [-a, 0, a],
        [a, a, 0],
        [-a, a, 0],
        [b, c, c],
        [-b, c, c],
        [b, -c, c],
        [-b, -c, c],
        [c, b, c],
        [-c, b, c],
        [c, -b, c],
        [-c, -b, c],
        [c, c, b],
        [-c, c, b],
        [c, -c, b],
        [-c, -c, b],
    ]
    w1 = 0.02652142440932
    w2 = 0.0199301476312
    w3 = 0.0250712367487

    w = 2 .* [fill(w1, 3); fill(w2, 6); fill(w3, 12)] # Multiply by two since integration is over the half-sphere

    ℒinv(x) = x * (3 - x^2) / (1 - x^2)

    function W(λ⃗)
        F = diagm(λ⃗)
        @tullio t⃗[i] := F * r⃗[i]
        @tullio n⃗[i] := inv(F') * r⃗[i]
        @tullio λ̄[i] := norm(t⃗[i])
        @tullio ν̄[i] := norm(n⃗[i])
        @tullio λ := (λ̄[i]^p) * w[i]# |> Base.Fix2(^, (1 / p))
        λr = λ^(1 / p) / √N
        β = ℒinv(λr)
        @tullio ν := ν̄[i]^q * w[i]# |> Base.Fix2(^, 1 / q)
        return N * U * μ * ν^(1 / q) + N * μ * (λr * β + log(β / sinh(β)))
    end
end

function parameters(ψ::NonaffineMicroSphere)
    return (:μ, :N, :p, :U, :q)
end

"""
Affine Micro-Sphere [^1]

Parameters: μ, N, p, U, q

Model: See Paper

---
[^1]: > Miehe C, Göktepe S, Lulei F. A micro-macro approach to rubber-like materials—part I: the non-affine micro-sphere model of rubber elasticity. Journal of the Mechanics and Physics of Solids. 2004 Nov 1;52(11):2617-60.
"""
struct AffineMicroSphere end

function StrainEnergyDensityFunction(ψ::AffineMicroSphere, (; μ, N, p, U, q))
    a = √(2) / 2
    b = 0.836095596749
    c = 0.387907304067
    r⃗ = [
        [0, 0, 1],
        [0, 1, 0],
        [1, 0, 0],
        [0, a, a],
        [0, -a, a],
        [a, 0, a],
        [-a, 0, a],
        [a, a, 0],
        [-a, a, 0],
        [b, c, c],
        [-b, c, c],
        [b, -c, c],
        [-b, -c, c],
        [c, b, c],
        [-c, b, c],
        [c, -b, c],
        [-c, -b, c],
        [c, c, b],
        [-c, c, b],
        [c, -c, b],
        [-c, -c, b],
    ]
    w1 = 0.0265214244093
    w2 = 0.0199301476312
    w3 = 0.0250712367487

    w = 2 .* [fill(w1, 3); fill(w2, 6); fill(w3, 12)] # Multiply by two since integration is over the half-sphere

    ℒinv(x) = x * (3 - x^2) / (1 - x^2)

    function W(λ⃗)
        F = diagm(λ⃗)
        @tullio t⃗[i] := F * r⃗[i]
        @tullio n⃗[i] := inv(F') * r⃗[i]
        @tullio λ̄[i] := norm(t⃗[i])
        @tullio ν̄[i] := norm(n⃗[i])
        @tullio λ := (λ̄[i]) * w[i]# |> Base.Fix2(^, (1 / p))
        λr = λ^(1 / p) / √N
        β = ℒinv(λr)
        @tullio ν := ν̄[i]^q * w[i]# |> Base.Fix2(^, 1 / q)
        return N * U * μ * ν^(1 / q) + N * μ * (λr * β + log(β / sinh(β)))
    end
end

function parameters(ψ::AffineMicroSphere)
    return (:μ, :N, :p, :U, :q)
end

"""
Bootstrapped 8Chain Model

Parameters: μ, N

Model: ``W_8(\\frac{\\sum\\lambda}{\\sqrt{3N}}-\\frac{\\lambda_{chain}}{\\sqrt{N}})+W_{8}(\\frac{\\lambda_{chain}}{\\sqrt{N}})``

``W_8(x) = \\mu N (x \\mathcal{L}^{-1}(x) + \\log\\frac{\\mathcal{L}^{-1}(x)}{\\sinh\\mathcal{L}^{-1}(x)})``

``\\lambda_{chain} = \\sqrt{\\frac{I_1}{3}}``
"""
struct Bootstrapped8Chain end

function StrainEnergyDensityFunction(ψ::Bootstrapped8Chain, (; μ, N))
    ℒinv(x) = x * (3 - 1.0651 * x^2 - 0.962245 * x^4 + 1.47353 * x^6 - 0.48953 * x^8) / (1 - x) / (1 + 1.01524 * x)
    W8(x) = μ * N * (x * ℒinv(x) + log(ℒinv(x) / sinh(ℒinv(x))))
    function W(λ⃗)
        λchain = √(I₁(λ⃗) / 3)
        W8(sum(λ⃗) / √(3N) - λchain / √(N)) + W8(λchain / √(N))
    end
end

function parameters(ψ::Bootstrapped8Chain)
    return (:μ, :N)
end

"""
Davidson - Goulbourne

Parameters: Gc, Ge, λmax

Model: ``\\frac{G_c}{6}I_1-G_c\\lambda_{max}\\log\\bigg(3\\lambda_{max}^2-I_1\\bigg)+G_e\\sum\\limits_{i=1}^{3}\\big(\\lambda_i+\\frac{1}{\\lambda_i}\\big)``
"""
struct DavidsonGoulbourne end

function StrainEnergyDensityFunction(ψ::DavidsonGoulbourne, (; Gc, Ge, λmax))
    W(λ⃗) = 1 / 6 * Gc * I₁(λ⃗) - Gc * λmax^2 * log(3λmax^2 - I₁(λ⃗)) + Ge * sum(λ⃗ .+ 1 ./ λ⃗)
end

function parameters(ψ::DavidsonGoulbourne)
    return (:Gc, :Ge, :λmax)
end

"""
Khiêm-Itskov Model [^1]

Parameters: μcκ, n, q, μt

Model: ``\\mu_c \\kappa n \\log\\bigg(\\frac{\\sin(\\frac{\\pi}{\\sqrt{n}})(\\frac{I_1}{3})^{\\frac{q}{2}}}{\\sin(\\frac{\\pi}{\\sqrt{n}}(\\frac{I_1}{3})^{\\frac{q}{2}}}\\bigg)+\\mu_t\\big[\\frac{I_2}{3}^{1/2} - 1 \\big]``

[^1]: > Khiêm VN, Itskov M. Analytical network-averaging of the tube model:: Rubber elasticity. Journal of the Mechanics and Physics of Solids. 2016 Oct 1;95:254-69.
"""
struct KhiemItskov end

function StrainEnergyDensityFunction(ψ::KhiemItskov, (; μcκ, n, q, μt))
    W(λ⃗) = μcκ * n * log((sin(π / sqrt(n)) * (I₁(λ⃗) / 3)^(q / 2)) / (sin(π / sqrt(n) * (I₁(λ⃗) / 3)^(q / 2)))) + μt * ((I₂(λ⃗) / 3)^(1 / 2) - 1)
end

function StrainEnergyDensityFunction(ψ::KhiemItskov, (; μcκ, n, q, μt), I::InvariantForm)
    W(I⃗) = μcκ * n * log((sin(π / sqrt(n)) * (I⃗[1] / 3)^(q / 2)) / (sin(π / sqrt(n) * (I⃗[1] / 3)^(q / 2)))) + μt * ((I⃗[2] / 3)^(1 / 2) - 1)
end

function parameters(ψ::KhiemItskov)
    return (:μcκ, :n, :q, :μt)
end

"""
General Constitutive Model

Parameters: Gc, Ge, N

Model: ``G_c N \\log\\bigg(\\frac{3N+\\frac{1}{2}I_1}{3N-I_1}\\bigg)+G_e\\sum\\limits_{i=1}^{3}\\frac{1}{\\lambda_I}``
"""
struct GeneralConstitutiveModel end

function StrainEnergyDensityFunction(ψ::GeneralConstitutiveModel, (; Gc, Ge, N))
    W(λ⃗) = Gc * N * log((3N + 0.5 * I₁(λ⃗)) / (3N - I₁(λ⃗))) + Ge * sum(λ⃗ .^ (-1))
end

function parameters(ψ::GeneralConstitutiveModel)
    return (:Gc, :Ge, :N)
end

"""
Full Network - Wu Geisson

Parameters: μ, N, ρ

Model: ``(1-\\rho)W_{3Chain}+\\rho W_{8chain}``
"""
struct FullNetwork end

function StrainEnergyDensityFunction(ψ::FullNetwork, (; μ, N, ρ))
    W3 = StrainEnergyDensityFunction(ThreeChainModel(), (μ=μ, N=N))
    W8 = StrainEnergyDensityFunction(ArrudaBoyce(), (μ=μ, N=N))
    W(λ⃗) = (1 - ρ) * W3(λ⃗) + ρ * W8(λ⃗)
end

function parameters(ψ::FullNetwork)
    return (:μ, :N, :ρ)
end

"""
Zuniga - Beatty

Parameters: μ, N₃, N₈

Model: ``\\sqrt{\\frac{N_3+N_8}{2N_3}}W_{3Chain}+\\sqrt{\\frac{I_1}{3N_8}}W_{8Chain}``
"""
struct ZunigaBeatty end

function StrainEnergyDensityFunction(ψ::ZunigaBeatty, (; μ, N₃, N₈))
    ΛL = √((N₃ + N₈) / 2)
    Λch = 1 / √(3) * √(I₁(λ⃗))
    ρ₃ = ΛL / √(N₃)
    ρ₈ = Λch / √(N₈)
    W3 = StrainEnergyDensityFunction(ThreeChainModel(), (μ=μ, N=N₃))
    W8 = StrainEnergyDensityFunction(ArrudaBoyce(), (μ=μ, N₈=N₈))
    W(λ⃗) = ρ₃ * W3(λ⃗) + ρ₈ * W8(λ⃗)
end

function parameters(ψ::ZunigaBeatty)
    return (:μ, :N₃, :N₈)
end

"""
Lim

Parameters: μ₁, μ₂, N, Î₁

Model: ``(1-f(\\frac{I_1-3}{\\hat{I_1}-3}))W_{NeoHookean}(μ₁)+fW_{ArrudaBoyce}(μ₂, N)``
"""
struct Lim end

function StrainEnergyDensityFunction(ψ::Lim, (; μ₁, μ₂, N, Î₁))
    Wg = StrainEnergyDensityFunction(NeoHookean(), (μ = μ₁))
    W8 = StrainEnergyDensityFunction(ArrudaBoyce(), (μ=μ₂, N=N))
    f(x) = x^3 * (10 - 15x + 6x^2)
    function W(λ⃗)
        ζ = (I₁(λ⃗) - 3) / (Î₁ - 3)
        (1 - f(ζ)) * Wg(λ⃗) + f(ζ) * W8(λ⃗)
    end
end

function StrainEnergyDensityFunction(ψ::Lim, (; μ₁, μ₂, N, Î₁), I::InvariantForm)
    Wg = StrainEnergyDensityFunction(NeoHookean(), (μ = μ₁), I)
    W8 = StrainEnergyDensityFunction(ArrudaBoyce(), (μ=μ₂, N=N), I)
    f(x) = x^3 * (10 - 15x + 6x^2)
    function W(I⃗)
        ζ = (I⃗[1] - 3) / (Î₁ - 3)
        (1 - f(ζ)) * Wg(I⃗) + f(ζ) * W8(I⃗)
    end
end

function parameters(ψ::Lim)
    return (:μ₁, :μ₂, :N, :Î₁)
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
struct BechirChevalier end

function StrainEnergyDensityFunction(ψ::BechirChevalier, (; μ₀, η, ρ, N₃, N₈))
    μf = ρ * √(I₁ / 3 / N₈)
    W3 = StrainEnergyDensityFunction(ThreeChainModel(), (μ=μf, N=N₃))
    function W(λ⃗)
        α = maximum(λ⃗)
        μc = (1 - η * α / √(N₃)) * μ₀
        W8 = StrainEnergyDensityFunction(ArrudaBoyce(), (μ=μc / 3, N=N₈))
        W3(λ⃗) + W8(λ⃗)
    end
end

function parameters(ψ::BechirChevalier)
    return (:μ₀, :η, :ρ, :N₃, :N₈)
end
