# # Available Models
export GeneralMooneyRivlin, GeneralDarijaniNaghdabadi, GeneralBeda, MooneyRivlin, NeoHookean, Gent, Biderman, Isihara, JamesGreenSimpson, Lion, Yeoh, HauptSedlan, HartmannNeff, HainesWilson, Carroll, BahremanDarijani, Zhao, Knowles, Swanson, YamashitaKawabata, DavisDeThomas, Gregory, ModifiedGregory, Beda, Amin, LopezPamies, GenYeoh, VerondaWestmann, FungDemiray, Vito, ModifiedYeoh, MansouriDarijani, GentThomas, HossMarczakI, HossMarczakII, ExpLn, VanDerWaals, TakamizawaHayashi, YeohFleming, PucciSaccomandi, HorganSaccomandi, Beatty, HorganMurphy, ArrudaBoyce, Ogden, EdwardVilgis, NonaffineTube, Tube, MCC, Bechir4Term, ConstrainedJunction, ContinuumHybrid, ArmanNarooei, PengLandel, ValanisLandel, Attard, Shariff, ThreeChainModel, ModifiedFloryErman, ABGI, BechirChevalier, Bootstrapped8Chain, DavidsonGoulbourne, ExtendedTubeModel, FullNetwork, HartSmith, GeneralConstitutiveModel, Lim, NonaffineMicroSphere, AffineMicroSphere, KhiemItskov, ZunigaBeatty, ChevalierMarco, Alexander, GornetDesmorat, LambertDianiRey

"""
General Mooney Rivlin[^1]

Parameters: [C]

Model:
``\\sum\\limits_{i,j = 0}^{N,M} C_{i,j}(I_1-3)^i(I_2-3)^j``

[^1]: > Mooney M. A theory of large elastic deformation. Journal of applied physics. 1940 Sep;11(9):582-92.
"""
struct GeneralMooneyRivlin <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::GeneralMooneyRivlin, λ⃗::AbstractVector, (; C))
    I1 = I₁(λ⃗)
    I2 = I₂(λ⃗)
    @tullio W := C[j, i] * (I1 - 3)^(i - 1) * (I2 - 3)^(j - 1)
end

function ContinuumModels.StrainEnergyDensity(ψ::GeneralMooneyRivlin, I⃗, (; C), I::InvariantForm)
    @tullio W := C[j, i] * (I⃗[1] - 3)^(i - 1) * (I⃗[2] - 3)^(j - 1)
    return W
end

parameters(ψ::GeneralMooneyRivlin) = (:C,)

citation(ψ::GeneralMooneyRivlin) = get_citation("mooney1940theory")

"""
General Darijani Naghdabadi [^1]

Parameters: A⃗, B⃗, m⃗, n⃗

Model: ``\\sum\\limits_{i = 1}{3}\\sum\\limits_{j=0}^{N} A_j (\\lambda_i^{m_j}-1) + B_j(\\lambda_i^{-n_j}-1)``

[^1]: > Bahreman M, Darijani H. New polynomial strain energy function; application to rubbery circular cylinders under finite extension and torsion. Journal of Applied Polymer Science. 2015 Apr 5;132(13).
"""
struct GeneralDarijaniNaghdabadi <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::GeneralDarijaniNaghdabadi, λ⃗::AbstractVector, (; A⃗, B⃗, m⃗, n⃗))
    @assert length(A⃗) == length(m⃗) "Length of A⃗ ≠ length of m⃗"
    @assert length(B⃗) == length(n⃗) "Length of B⃗ ≠ length of n⃗"
    sum(i -> sum(A⃗ .* (λ⃗[i] .^ m⃗ .- 1)) + sum(B⃗ .* (λ⃗[i] .^ (-1 .* n⃗) .- 1)), 1:3)
    # sum(i -> sum(A⃗ .* (λ⃗[i] .^ m⃗ .- 1)) + sum(B⃗ .* (λ⃗[i] .^ (-1 .* n⃗) .- 1)), [1:3])
end

function parameters(ψ::GeneralDarijaniNaghdabadi)
    return (:A⃗, :B⃗, :m⃗, :n⃗)
end

citation(ψ::GeneralDarijaniNaghdabadi) = get_citation("bahreman2015new")

"""
General Beda [^1]

Parameters: C, K, α, β

Model: ``\\sum\\limits_{i = 1}^{N}\\frac{C_i}{\\alpha_i}(I_1-3)^{\\alpha_i} + \\sum\\limits_{j=1}^{M}\\frac{K_j}{\\beta_j}(I_2-3)^{\\beta_j}``

[^1]: > Beda T. Reconciling the fundamental phenomenological expression of the strain energy of rubber with established experimental facts. Journal of Polymer Science Part B: Polymer Physics. 2005 Jan 15;43(2):125-34.
"""
struct GeneralBeda <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::GeneralBeda, λ⃗::AbstractVector, (; C, K, α, β))
    @assert length(C) == length(α) "Vector C and Vector α are not the same length"
    @assert length(K) == length(β) "Vector K and Vector β are not the same length"
    W1 = C ./ α .* (I₁(λ⃗) - 3) .^ α |> sum
    W2 = K ./ β .* (I₂(λ⃗) - 3) .^ β |> sum
    return W1 + W2
end

function ContinuumModels.StrainEnergyDensity(ψ::GeneralBeda, I⃗, (; C, K, α, β), I::InvariantForm)
    @assert length(C) == length(α) "Vector C and Vector α are not the same length"
    @assert length(K) == length(β) "Vector K and Vector β are not the same length"
    W1 = C ./ α .* (I⃗[1] - 3) .^ α |> sum
    W2 = K ./ β .* (I⃗[2] - 3) .^ β |> sum
    return W1 + W2
end

parameters(ψ::GeneralBeda) = (:C, :K, :α, :β)

citation(ψ::GeneralBeda) = get_citation("beda2005reconciling")

"""
Mooney Rivlin Model [^1]

Parameters: C01, C10

Model: ``C_{10}(I_1-3)+C_{01}(I_2-3)``

[^1]: > Mooney M. A theory of large elastic deformation. Journal of applied physics. 1940 Sep;11(9):582-92.
"""
struct MooneyRivlin <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::MooneyRivlin, λ⃗::AbstractVector, (; C10, C01))
    ContinuumModels.StrainEnergyDensity(
        GeneralMooneyRivlin(),
        λ⃗,
        (C=[
            0 C10
            C01 0
        ],
        )
    )
end

function ContinuumModels.StrainEnergyDensity(ψ::MooneyRivlin, I⃗, (; C10, C01), I::InvariantForm)
    ContinuumModels.StrainEnergyDensity(
        GeneralMooneyRivlin(),
        I⃗,
        (C=[
            0 C10
            C01 0
        ],
        ),
        I
    )
end

parameters(ψ::MooneyRivlin) = (:C10, :C01)

citation(ψ::MooneyRivlin) = get_citation("mooney1940theory")

"""
NeoHookean [^1]

Parameters: μ

Model: ``\\frac{\\mu}{2}(I_1-3)``

[^1]: > Treloar LR. The elasticity of a network of long-chain molecules—II. Transactions of the Faraday Society. 1943;39:241-6.
"""
struct NeoHookean <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::NeoHookean, λ⃗::AbstractVector, (; μ))
    μ / 2 * (I₁(λ⃗) - 3)
end

function ContinuumModels.StrainEnergyDensity(ψ::NeoHookean, I⃗, (; μ), I::InvariantForm)
    μ / 2 * (I⃗[1] - 3)
end

parameters(ψ::NeoHookean) = (:μ,)

citation(ψ::NeoHookean) = get_citation("treloar1943elasticity")

"""
Isihara [^1]

Parameters: C10, C20, C01

Model: ``\\sum\\limits_{i,j=0}^{2, 1}C_{i,j}(I_1-3)^i(I_2-3)^j``

[^1]: > Isihara A, Hashitsume N, Tatibana M. Statistical theory of rubber‐like elasticity. IV.(two‐dimensional stretching). The Journal of Chemical Physics. 1951 Dec;19(12):1508-12.
"""
struct Isihara <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Isihara, λ⃗::AbstractVector, (; C10, C20, C01))
    ContinuumModels.StrainEnergyDensity(
        GeneralMooneyRivlin(),
        λ⃗,
        (C=[
            0 C10 C20
            C01 0 0
        ],
        )
    )
end

function ContinuumModels.StrainEnergyDensity(ψ::Isihara, I⃗, (; C10, C20, C01), I::InvariantForm)
    ContinuumModels.StrainEnergyDensity(
        GeneralMooneyRivlin(),
        I⃗,
        (C=[
            0 C10 C20
            C01 0 0
        ],
        ),
        I
    )
end

parameters(ψ::Isihara) = (:C10, :C20, :C01)

citation(ψ::Isihara) = get_citation("isihara1951statistical")

"""
Biderman [^1]

Parameters: C10, C01, C20, C30

Model: ``\\sum\\limits_{i,j=0}^{3, 1}C_{i,j}(I_1-3)^i(I_2-3)^j``

[^1]: > Biderman VL. Calculation of rubber parts. Rascheti na prochnost. 1958;40.
"""
struct Biderman <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Biderman, λ⃗::AbstractVector, (; C10, C01, C20, C30))
    ContinuumModels.StrainEnergyDensity(
        GeneralMooneyRivlin(),
        λ⃗,
        (C=[
            0 C10 C20 C30
            C01 0 0 0
        ],
        )
    )
end

function ContinuumModels.StrainEnergyDensity(ψ::Biderman, I⃗, (; C10, C01, C20, C30), I::InvariantForm)
    ContinuumModels.StrainEnergyDensity(
        GeneralMooneyRivlin(),
        I⃗,
        (C=[
            0 C10 C20 C30
            C01 0 0 0
        ],
        ),
        I
    )
end

parameters(ψ::Biderman) = (:C10, :C01, :C20, :C30)

citation(ψ::Biderman) = get_citation("biderman1958calculation")

"""
James-Green-Simpson [^1]

Parameters: C10, C01, C11, C20, C30

Model: ``\\sum\\limits_{i,j=0}^{3, 1}C_{i,j}(I_1-3)^i(I_2-3)^j``

[^1]: > James AG, Green A, Simpson GM. Strain energy functions of rubber. I. Characterization of gum vulcanizates. Journal of Applied Polymer Science. 1975 Jul;19(7):2033-58.
"""
struct JamesGreenSimpson <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::JamesGreenSimpson, λ⃗::AbstractVector, (; C10, C01, C11, C20, C30))
    ContinuumModels.StrainEnergyDensity(
        GeneralMooneyRivlin(),
        λ⃗,
        (C=[
            0 C10 C20 C30
            C01 0 0 0
        ],
        )
    )
end

function ContinuumModels.StrainEnergyDensity(ψ::JamesGreenSimpson, I⃗, (; C10, C01, C11, C20, C30), I::InvariantForm)
    ContinuumModels.StrainEnergyDensity(
        GeneralMooneyRivlin(),
        I⃗,
        (C=[
            0 C10 C20 C30
            C01 0 0 0
        ],
        ),
        I
    )
end

parameters(ψ::JamesGreenSimpson) = (:C10, :C01, :C11, :C20, :C30)

citation(ψ::JamesGreenSimpson) = get_citation("james1975strain")

"""
Haines-Wilson [^1]

Parameters: C10, C01, C11, C02, C20, C30

Model: ``\\sum\\limits_{i,j=0}^{3, 2}C_{i,j}(I_1-3)^i(I_2-3)^j``

[^1]: > Haines DW, Wilson WD. Strain-energy density function for rubberlike materials. Journal of the Mechanics and Physics of Solids. 1979 Aug 1;27(4):345-60.
"""
struct HainesWilson <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::HainesWilson, λ⃗::AbstractVector, (; C10, C01, C11, C02, C20, C30))
    ContinuumModels.StrainEnergyDensity(
        GeneralMooneyRivlin(),
        λ⃗,
        (C=[
            0 C10 C20 C30
            C01 C11 0 0
            C02 0 0 0
        ],
        )
    )
end

function ContinuumModels.StrainEnergyDensity(ψ::HainesWilson, I⃗, (; C10, C01, C11, C02, C20, C30), I::InvariantForm)
    ContinuumModels.StrainEnergyDensity(
        GeneralMooneyRivlin(),
        I⃗,
        (C=[
            0 C10 C20 C30
            C01 C11 0 0
            C02 0 0 0
        ],
        ),
        I
    )
end

parameters(ψ::HainesWilson) = (:C10, :C01, :C11, :C02, :C20, :C30)

citation(ψ::HainesWilson) = get_citation("haines1979strain")

"""
Yeoh [^1]

Parameters: C10, C20, C30

Model: ``\\sum\\limits_{i,j=0}^{3, 0}C_{i,j}(I_1-3)^i(I_2-3)^j``

[^1]: > Yeoh OH. Characterization of elastic properties of carbon-black-filled rubber vulcanizates. Rubber chemistry and technology. 1990 Nov;63(5):792-805.
"""
struct Yeoh <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Yeoh, λ⃗::AbstractVector, (; C10, C20, C30))
    ContinuumModels.StrainEnergyDensity(
        GeneralMooneyRivlin(),
        λ⃗,
        (C=[0 C10 C20 C30],)
    )
end

function ContinuumModels.StrainEnergyDensity(ψ::Yeoh, I⃗, (; C10, C20, C30), I::InvariantForm)
    ContinuumModels.StrainEnergyDensity(
        GeneralMooneyRivlin(),
        I⃗,
        (C=[0 C10 C20 C30],),
        I
    )
end

parameters(ψ::Yeoh) = (:C10, :C20, :C30)

citation(ψ::Yeoh) = get_citation("yeoh1990characterization")

"""
Lion [^1]

Parameters: C10, C01, C50

Model: ``\\sum\\limits_{i,j=0}^{5,1}C_{i,j}(I_1-3)^i(I_2-3)^j``

[^1]: > Lion A. On the large deformation behaviour of reinforced rubber at different temperatures. Journal of the Mechanics and Physics of Solids. 1997 Nov 1;45(11-12):1805-34.
"""
struct Lion <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Lion, λ⃗::AbstractVector, (; C10, C01, C50))
    ContinuumModels.StrainEnergyDensity(
        GeneralMooneyRivlin(),
        λ⃗,
        (C=[
            0 C10 0 0 0 C50
            C01 0 0 0 0 0
        ],)
    )
end

function ContinuumModels.StrainEnergyDensity(ψ::Lion, I⃗, (; C10, C01, C50), I::InvariantForm)
    ContinuumModels.StrainEnergyDensity(
        GeneralMooneyRivlin(),
        I⃗,
        (C=[
            0 C10 0 0 0 C50
            C01 0 0 0 0 0
        ],),
        I
    )
end

parameters(ψ::Lion) = (:C10, :C01, :C50)

citation(ψ::Lion) = get_citation("lion1997large")


"""
Haupt Sedlan [^1]

Parameters: C10, C01, C11, C02, C30

Model:
``\\sum\\limits_{i,j=0}^{3, 2}C_{i,j}(I_1-3)^i(I_2-3)^j``

[^1]: > Haupt P, Sedlan K. Viscoplasticity of elastomeric materials: experimental facts and constitutive modelling. Archive of Applied Mechanics. 2001 Mar;71(2):89-109.
"""
struct HauptSedlan <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::HauptSedlan, λ⃗::AbstractVector, (; C10, C01, C11, C02, C30))
    ContinuumModels.StrainEnergyDensity(
        GeneralMooneyRivlin(),
        λ⃗,
        (C=[
            0 C10 0 C30
            C01 C11 0 0
            C02 0 0 0
        ],)
    )
end

function ContinuumModels.StrainEnergyDensity(ψ::HauptSedlan, I⃗, (; C10, C01, C11, C02, C30), I::InvariantForm)
    ContinuumModels.StrainEnergyDensity(
        GeneralMooneyRivlin(),
        I⃗,
        (C=[
            0 C10 0 C30
            C01 C11 0 0
            C02 0 0 0
        ],),
        I
    )
end

parameters(ψ::HauptSedlan) = (:C10, :C01, :C11, :C02, :C30)

citation(ψ::HauptSedlan) = get_citation("haupt2001viscoplasticity")

"""
Hartmann-Neff [^1]

Parameters: α, Ci0, C0j

Model: ``\\sum\\limits_{i,j=0}^{M,N}C_{i,0}(I_1-3)^i -3\\sqrt{3}^j+\\alpha(I_1-3)``

[^1]: > Hartmann S, Neff P. Polyconvexity of generalized polynomial-type hyperelastic strain energy functions for near-incompressibility. International journal of solids and structures. 2003 Jun 1;40(11):2767-91.
"""
struct HartmannNeff <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::HartmannNeff, λ⃗::AbstractVector, (; α, Ci⃗0, C0j⃗))
    @tullio W1 := Ci⃗0[i] * (I₁(λ⃗) - 3)^i
    @tullio W2 := C0j⃗[j] * (I₂(λ⃗)^(3 / 2) - 3sqrt(3))^j
    return W1 + W2 + α * (I₁(λ⃗)^3 - 3^3)
end

function ContinuumModels.StrainEnergyDensity(ψ::HartmannNeff, I⃗, (; α, Ci⃗0, C0j⃗), I::InvariantForm)
    @tullio W1 := Ci⃗0[i] * (I⃗[1] - 3)^i
    @tullio W2 := C0j⃗[j] * (I⃗[2]^(3 / 2) - 3sqrt(3))^j
    return W1 + W2 + α * (I⃗[1]^3 - 3^3)
end

parameters(ψ::HartmannNeff) = (:α, :Ci⃗0, :C0j⃗)

citation(ψ::HartmannNeff) = get_citation("carroll2011strain")

"""
Carroll [^1]

Parameters: A, B, C

Model: ``AI_1+BI_1^4+C\\sqrt{I_2}``

[^1]: > Carroll M. A strain energy function for vulcanized rubbers. Journal of Elasticity. 2011 Apr;103(2):173-87.
"""
struct Carroll <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Carroll, λ⃗::AbstractVector, (; A, B, C))
    A * I₁(λ⃗) + B * I₁(λ⃗)^4 + C * I₂(λ⃗)^(1 / 2)
end

function ContinuumModels.StrainEnergyDensity(ψ::Carroll, I⃗, (; A, B, C), I::InvariantForm)
    A * I⃗[1] + B * I⃗[1]^4 + C * I⃗[2]^(1 / 2)
end

parameters(ψ::Carroll) = (:A, :B, :C)

citation(ψ::Carroll) = get_citation("carroll2011strain")

"""
Bahreman Darijani [^1]

Parameters: A2, B2, A4, A6

Model:
``\\sum\\limits_{i = 1}{3}\\sum\\limits_{j=0}^{N} A_j (\\lambda_i^{m_j}-1) + B_j(\\lambda_i^{-n_j}-1)``

[^1]: > Bahreman M, Darijani H. New polynomial strain energy function; application to rubbery circular cylinders under finite extension and torsion. Journal of Applied Polymer Science. 2015 Apr 5;132(13).
"""
struct BahremanDarijani <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::BahremanDarijani, λ⃗::AbstractVector, (; A2, B2, A4, A6))
    ContinuumModels.StrainEnergyDensity(
        GeneralDarijaniNaghdabadi(),
        λ⃗,
        (
            A⃗=[0, A2, 0, A4, 0, A6],
            B⃗=[0, B2],
            m⃗=[0, 2, 0, 4, 0, 6],
            n⃗=[0, 2])
    )
end

parameters(ψ::BahremanDarijani) = (:A2, :B2, :A4, :A6)

citation(ψ::BahremanDarijani) = get_citation("bahreman2015new")

"""
Zhao [^1]

Parameters: C₋₁¹,, C₁¹, C₂¹, C₂²

Model: ``C_{-1}^1*(I_2-3)+C_{1}^{1}(I_1-3)+C_{2}^{1}(I_1^2-2I_2-3)+C_{2}^{2}(I_1^2-2I_2-3)^2``

[^1]: > Zhao Z, Mu X, Du F. Modeling and verification of a new hyperelastic model for rubber-like materials. Mathematical Problems in Engineering. 2019 May 2;2019.
"""
struct Zhao <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Zhao, λ⃗::AbstractVector, (; C₋₁¹, C₁¹, C₂¹, C₂²))
    C₋₁¹ * (I₂(λ⃗) - 3) + C₁¹ * (I₁(λ⃗) - 3) + C₂¹ * (I₁(λ⃗)^2 - 2I₂(λ⃗) - 3) + C₂² * (I₁(λ⃗)^2 - 2I₂(λ⃗) - 3)^2
end

function ContinuumModels.StrainEnergyDensity(ψ::Zhao, (; C₋₁¹, C₁¹, C₂¹, C₂²), I::InvariantForm)
    C₋₁¹ * (I⃗[2] - 3) + C₁¹ * (I⃗[1] - 3) + C₂¹ * (I⃗[1]^2 - 2I⃗[2] - 3) + C₂² * (I⃗[1]^2 - 2I⃗[2] - 3)^2
end

parameters(ψ::Zhao) = (:C₋₁¹, :C₁¹, :C₂¹, :C₂²)

citation(ψ::Zhao) = get_citation("zhao2019modeling")

"""
Knowles [^1]

Parameters: μ, b, n

Model: ``\\frac{\\mu}{2b}((1+\\frac{b}{n}(I_1-3))^n-1)``

[^1]: > Knowles JK. The finite anti-plane shear field near the tip of a crack for a class of incompressible elastic solids. International Journal of Fracture. 1977 Oct;13(5):611-39.
"""
struct Knowles <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Knowles, λ⃗::AbstractVector, (; μ, b, n))
    μ / (2b) * ((1 + (b / n) * (I₁(λ⃗) - 3))^n - 1)
end

function ContinuumModels.StrainEnergyDensity(ψ::Knowles, I⃗, (; μ, b, n), I::InvariantForm)
    μ / (2b) * ((1 + (b / n) * (I⃗[1] - 3))^n - 1)
end


parameters(ψ::Knowles) = (:μ, :b, :n)

function parameter_bounds(ψ::Knowles, data::AbstractHyperelasticData)
    lb = (μ=-Inf, b=0, n=0)
    ub = nothing
    return (lb=lb, ub=ub)
end

citation(ψ::Knowles) = get_citation("knowles1977finite")

"""
Swanson [^1]

Parameters: A⃗, α⃗, B⃗, β⃗

Model: ``\\sum\\limits_{i=1}^{N} \\frac{3}{2}(\\frac{A_i}{1+\\alpha_i}(\\frac{I_1}{3})^{1+\\alpha_i}+\\frac{B_i}{1+\\beta_i}(\\frac{I_2}{3})^{1+\\beta_i}``

[^1]: > Swanson SR. A constitutive model for high elongation elastic materials.
"""
struct Swanson <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Swanson, λ⃗::AbstractVector, (; A⃗, α⃗, B⃗, β⃗))
    @assert length(A⃗) == length(α⃗) == length(B⃗) == length(β⃗) "The vectors are not the same length"
    @tullio _ := 3 / 2 * (A⃗[i] / (1 + α⃗[i]) * (I₁(λ⃗) / 3)^(1 + α⃗[i]) + B⃗[i] / (1 + β⃗[i]) * (I₂(λ⃗) / 3)^(1 + β⃗[i]))
end

function ContinuumModels.StrainEnergyDensity(ψ::Swanson, I⃗, (; A⃗, α⃗, B⃗, β⃗), I::InvariantForm)
    @assert length(A⃗) == length(α⃗) == length(B⃗) == length(β⃗) "The vectors are not the same length"
    @tullio _ := 3 / 2 * (A⃗[i] / (1 + α⃗[i]) * (I⃗[1] / 3)^(1 + α⃗[i]) + B⃗[i] / (1 + β⃗[i]) * (I⃗[2] / 3)^(1 + β⃗[i]))
end

parameters(ψ::Swanson) = (:A⃗, :α⃗, :B⃗, :β⃗)

citation(ψ::Swanson) = get_citation("swanson1985constitutive")

"""
Yamashita-Kawabata [^1]

Parameters: C1, C2, C3, N

Model: ``C_1(I_1-3)+C_2(I_2-3)+\\frac{C_3}{N+1}(I_1-3)^{N+1}``

[^1]: > Yamashita Y, Kawabata S. Approximated form of the strain energy-density function of carbon-black filled rubbers for industrial applications. Nippon Gomu Kyokaishi(Journal of the Society of Rubber Industry, Japan)(Japan). 1992;65(9):517-28.
"""
struct YamashitaKawabata <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::YamashitaKawabata, λ⃗::AbstractVector, (; C1, C2, C3, N))
    C1 * (I₁(λ⃗) - 3) + C2 * (I₂(λ⃗) - 3) + C3 / (N + 1) * (I₁(λ⃗) - 3)^(N + 1)
end

function ContinuumModels.StrainEnergyDensity(ψ::YamashitaKawabata, I⃗, (; C1, C2, C3, N), I::InvariantForm)
    1 * (I⃗[1] - 3) + C2 * (I⃗[2] - 3) + C3 / (N + 1) * (I⃗[1] - 3)^(N + 1)
end

parameters(ψ::YamashitaKawabata) = (:C1, :C2, :C3, :N)

citation(ψ::YamashitaKawabata) = get_citation("yamashita1992approximated")

"""
Davis-DeThomas [^1]

Parameters: A, n, C, k

Model: ``\\frac{A}{2(1-\\frac{n}{2})}(I_1-3+C^2)^{1-\\frac{n}{2}}+k(I_1-3)^2``

[^1]: > Davies CK, De DK, Thomas AG. Characterization of the behavior of rubber for engineering design purposes. 1. Stress-strain relations. Rubber chemistry and technology. 1994 Sep;67(4):716-28.
"""
struct DavisDeThomas <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::DavisDeThomas, λ⃗::AbstractVector, (; A, n, C, k))
    A / (2 * (1 - n / 2)) * (I₁(λ⃗) - 3 + C^2)^(1 - n / 2) + k * (I₁(λ⃗) - 3)^2
end

function ContinuumModels.StrainEnergyDensity(ψ::DavisDeThomas, I⃗, (; A, n, C, k), I::InvariantForm)
    A / (2 * (1 - n / 2)) * (I⃗[1] - 3 + C^2)^(1 - n / 2) + k * (I⃗[1] - 3)^2
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
struct Gregory <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Gregory, λ⃗::AbstractVector, (; A, B, C, m, n))
    A / (2 - n) * (I₁(λ⃗) - 3 + C^2)^(1 - n / 2) + B / (2 + m) * (I₁(λ⃗) - 3 + C^2)^(1 + m / 2)
end

function ContinuumModels.StrainEnergyDensity(ψ::Gregory, I⃗, (; A, B, C, m, n), I::InvariantForm)
    A / (2 - n) * (I⃗[1] - 3 + C^2)^(1 - n / 2) + B / (2 + m) * (I⃗[1] - 3 + C^2)^(1 + m / 2)
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
struct ModifiedGregory <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::ModifiedGregory, λ⃗::AbstractVector, (; A, α, M, B, β, N))
    A / (1 + α) * (I₁(λ⃗) - 3 + M^2)^(1 + α) + B / (1 + β) * (I₁(λ⃗) - 3 + N^2)^(1 + β)
end

function ContinuumModels.StrainEnergyDensity(ψ::ModifiedGregory, I⃗, (; A, α, M, B, β, N), I::InvariantForm)
    A / (1 + α) * (I⃗[1] - 3 + M^2)^(1 + α) + B / (1 + β) * (I⃗[1] - 3 + N^2)^(1 + β)
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
struct Beda <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Beda, λ⃗::AbstractVector, (; C1, C2, C3, K1, α, β, ζ))
    ContinuumModels.StrainEnergyDensity(
        GeneralBeda(),
        λ⃗,
        (
            C=[C1, C2, C3],
            K=[K1],
            α=[α, 1, ζ],
            β=[β]
        )
    )
end

function ContinuumModels.StrainEnergyDensity(ψ::Beda, I⃗, (; C1, C2, C3, K1, α, β, ζ), I::InvariantForm)
    ContinuumModels.StrainEnergyDensity(
        GeneralBeda(),
        I⃗,
        (
            C=[C1, C2, C3],
            K=[K1],
            α=[α, 1, ζ],
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
struct Amin <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Amin, λ⃗::AbstractVector, (; C1, C2, C3, C4, N, M))
    C1 * (I₁(λ⃗) - 3) + C2 / (N + 1) * (I₁(λ⃗) - 3)^(N + 1) + C3 / (M + 1) * (I₁(λ⃗) - 3)^(M + 1) + C4 * (I₂(λ⃗) - 3)
end

function ContinuumModels.StrainEnergyDensity(ψ::Amin, I⃗, (; C1, C2, C3, C4, N, M), I::InvariantForm)
    C1 * (I⃗[1] - 3) + C2 / (N + 1) * (I⃗[1] - 3)^(N + 1) + C3 / (M + 1) * (I⃗[1] - 3)^(M + 1) + C4 * (I⃗[2] - 3)
end

function parameters(ψ::Amin)
    return (:C1, :C2, :C3, :C4, :N, :M)
end

"""
Lopez-Pamies [^1]

Parameters: α⃗, μ⃗

Model: ``\\frac{3^{1 - \\alpha_i}}{2\\alpha_i} \\mu_i (I_1^{\\alpha_i} - 3^{\\alpha_i})``

[^1]: > Lopez-Pamies O. A new I1-based hyperelastic model for rubber elastic materials. Comptes Rendus Mecanique. 2010 Jan 1;338(1):3-11.
"""
struct LopezPamies <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::LopezPamies, λ⃗::AbstractVector, (; α⃗, μ⃗))
    @assert length(α⃗) == length(μ⃗) "length of α⃗ is not equal to length of μ⃗"
    @tullio _ := (3^(1 - α⃗[i])) / (2α⃗[i]) * μ⃗[i] * (I₁(λ⃗)^(α⃗[i]) - 3^(α⃗[i]))
end

function ContinuumModels.StrainEnergyDensity(ψ::LopezPamies, I⃗, (; α⃗, μ⃗), I::InvariantForm)
    @assert length(α⃗) == length(μ⃗) "length of α⃗ is not equal to length of μ⃗"
    @tullio _ := (3^(1 - α⃗[i])) / (2α⃗[i]) * μ⃗[i] * (I⃗[1]^(α⃗[i]) - 3^(α⃗[i]))
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
struct GenYeoh <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::GenYeoh, λ⃗::AbstractVector, (; K1, K2, K3, m, p, q))
    K1 * (I₁(λ⃗) - 3)^m + K2 * (I₁(λ⃗) - 3)^p + K3 * (I₁(λ⃗) - 3)^q
end

function ContinuumModels.StrainEnergyDensity(ψ::GenYeoh, I⃗, (; K1, K2, K3, m, p, q), I::InvariantForm)
    K1 * (I⃗[1] - 3)^m + K2 * (I⃗[1] - 3)^p + K3 * (I⃗[1] - 3)^q
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
struct HartSmith <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::HartSmith, λ⃗::AbstractVector, (; G, k₁, k₂))
    G * exp(-9k₁ + k₁ * I₁(λ⃗)) / k₁ + G * k₂ * log(I₂(λ⃗))
end

function ContinuumModels.StrainEnergyDensity(ψ::HartSmith, I⃗, (; G, k₁, k₂), I::InvariantForm)
    G * exp(-9k₁ + k₁ * I⃗[1]) / k₁ + G * k₂ * log(I⃗[2])
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
struct VerondaWestmann <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::VerondaWestmann, λ⃗::AbstractVector, (; C1, C2, α))
    C1 * (exp(α * (I₁(λ⃗) - 3)) - 1) + C2 * (I₂(λ⃗) - 3)
end

function ContinuumModels.StrainEnergyDensity(ψ::VerondaWestmann, I⃗, (; C1, C2, α), I::InvariantForm)
    C1 * (exp(α * (I⃗[1] - 3)) - 1) + C2 * (I⃗[2] - 3)
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
struct FungDemiray <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::FungDemiray, λ⃗::AbstractVector, (; μ, b))
    μ / (2 * b) * (exp(b * (I₁(λ⃗) - 3)) - 1)
end

function ContinuumModels.StrainEnergyDensity(ψ::FungDemiray, I⃗, (; μ, b), I::InvariantForm)
    μ / (2 * b) * (exp(b * (I⃗[1] - 3)) - 1)
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
struct Vito <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Vito, λ⃗::AbstractVector, (; α, β, γ))
    α * (exp(β * (I₁(λ⃗) - 3) + γ * (I₂(λ⃗) - 3)) - 1)
end

function ContinuumModels.StrainEnergyDensity(ψ::Vito, I⃗, (; α, β, γ), I::InvariantForm)
    α * (exp(β * (I⃗[1] - 3) + γ * (I⃗[2] - 3)) - 1)
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
struct ModifiedYeoh <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::ModifiedYeoh, λ⃗::AbstractVector, (; C10, C20, C30, α, β))
    C10 * (I₁(λ⃗) - 3) + C20 * (I₁(λ⃗) - 3)^2 + C30 * (I₁(λ⃗) - 3)^3 + α / β * (1 - exp(-β * (I₁(λ⃗) - 3)))
end

function ContinuumModels.StrainEnergyDensity(ψ::ModifiedYeoh, I⃗, (; C10, C20, C30, α, β), I::InvariantForm)
    C10 * (I⃗[1] - 3) + C20 * (I⃗[1] - 3)^2 + C30 * (I⃗[1] - 3)^3 + α / β * (1 - exp(-β * (I⃗[1] - 3)))
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
struct ChevalierMarco <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::ChevalierMarco, λ⃗::AbstractVector, (; a⃗, b⃗))
    ∂W∂I1(I₁) = exp(sum(@tullio _ := a⃗[i] * (I₁ - 3)^(i - 1)))
    ∂W∂I2(I₂) = @tullio _ := b⃗[i] / I₂^(i - 1)
    quadgk(∂W∂I1, 3, I₁(λ⃗))[1] + quadgk(∂W∂I2, 3, I₂(λ⃗))[1]
end

function ContinuumModels.StrainEnergyDensity(ψ::ChevalierMarco, I⃗, (; a⃗, b⃗), I::InvariantForm)
    ∂W∂I1(I₁) = exp(sum(@tullio _ := a⃗[i] * (I₁ - 3)^(i - 1)))
    ∂W∂I2(I₂) = @tullio _ := b⃗[i] / I₂^(i - 1)
    quadgk(∂W∂I1, 3, I⃗[1])[1] + quadgk(∂W∂I2, 3, I⃗[2])[1]
end

function NominalStressFunction(ψ::ChevalierMarco, λ⃗::AbstractVector, (; a⃗, b⃗))
    ∂W∂I1(λ⃗) = exp(sum(@tullio _ := a⃗[i] * (I₁(λ⃗) - 3)^(i - 1)))
    ∂W∂I2(λ⃗) = @tullio _ := b⃗[i] / I₂(λ⃗)^(i - 1)
    𝐒 = 2 * (I(3) * ∂W∂I1 - diagm(λ⃗ .^ 2)^(-2) * ∂W∂I2)
    sᵢ = diag(𝐒)
    sᵢ = sᵢ .- sᵢ[3] .* λ⃗[3] / λ⃗[1]
    return sᵢ
end

function true_stress(ψ::ChevalierMarco, (; a⃗, b⃗))
    ∂W∂I1(λ⃗) = exp(sum(@tullio _ := a⃗[i] * (I₁(λ⃗) - 3)^(i - 1)))
    ∂W∂I2(λ⃗) = @tullio _ := b⃗[i] / I₂(λ⃗)^(i - 1)
    s(λ⃗) = NominalStressFunction(ψ, λ⃗, (a⃗=a⃗, b⃗=b⃗))
    σᵢ = map(λ⃗ᵢ -> λ⃗ᵢ .* s(λ⃗ᵢ), λ⃗)
    return σᵢ
end

function parameters(ψ::ChevalierMarco)
    return (:a⃗, :b⃗)
end

"""
Gornet - Desmorat [^1]

Parameters: h₁, h₂, h₃

Model: ``W = h_1\\int\\exp{h_3(I_1-3)^2}\\text{d}I_1+3h_2\\int\\frac{1}{\\sqrt{I_2}}\\text{d}I_2 = \\frac{h_1 \\sqrt{\\pi} \\text{erfi}(\\sqrt{h_3}(I_1-3)^2)}{2\\sqrt{h_3}}+6h_2\\sqrt{I_2}``

* Note: the differential form was original form and the closed form SEF was determine via symbolic integration in Mathematica.

[^1]: > Gornet L, Marckmann G, Desmorat R, Charrier P. A new isotropic hyperelastic strain energy function in terms of invariants and its derivation into a pseudo-elastic model for Mullins effect: application to finite element analysis. Constitutive Models for Rubbers VII. 2012:265-71.
"""
struct GornetDesmorat <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::GornetDesmorat, λ⃗::AbstractVector, (; h₁, h₂, h₃))
    h₁ * √π * erfi(√h₃ * (I₁(λ⃗) - 3)^2) / 2 / √h₃ + 6 * h₂ * √(I₂(λ⃗))
end

function ContinuumModels.StrainEnergyDensity(ψ::GornetDesmorat, I⃗, (; h₁, h₂, h₃), I::InvariantForm)
    h₁ * √π * erfi(√h₃ * (I⃗[1] - 3)^2) / 2 / √h₃ + 6 * h₂ * √(I⃗[2])
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
struct MansouriDarijani <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::MansouriDarijani, λ⃗::AbstractVector, (; A1, m1, B1, n1))
    A1 * (exp(m1 * (I₁(λ⃗) - 3)) - 1) + B1 * (exp(n1 * (I₂(λ⃗) - 3)) - 1)
end

function ContinuumModels.StrainEnergyDensity(ψ::MansouriDarijani, I⃗, (; A1, m1, B1, n1), I::InvariantForm)
    A1 * (exp(m1 * (I⃗[1] - 3)) - 1) + B1 * (exp(n1 * (I⃗[2] - 3)) - 1)
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
struct GentThomas <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::GentThomas, λ⃗::AbstractVector, (; C1, C2))
    C1 * (I₁(λ⃗) - 3) + C2 * log(I₂(λ⃗) / 3)
end

function ContinuumModels.StrainEnergyDensity(ψ::GentThomas, I⃗, (; C1, C2), I::InvariantForm)
    C1 * (I⃗[1] - 3) + C2 * log(I⃗[2] / 3)
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
struct Alexander <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Alexander, λ⃗::AbstractVector, (; C₁, C₂, C₃, k, γ))
    C₁ * √π * erfi(√k * (I₁(λ⃗) - 3)) / 2 / √k + C₂ * log((I₂(λ⃗) - 3 + γ) / γ) + C₃ * (I₂(λ⃗) - 3)
end

function ContinuumModels.StrainEnergyDensity(ψ::Alexander, I⃗, (; C₁, C₂, C₃, k, γ), I::InvariantForm)
    C₁ * √π * erfi(√k * (I⃗[1] - 3)) / 2 / √k + C₂ * log((I⃗[2] - 3 + γ) / γ) + C₃ * (I⃗[2] - 3)
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
struct LambertDianiRey <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::LambertDianiRey, λ⃗::AbstractVector, (; a⃗, b⃗))
    ∂W∂I₁(I₁) = exp(@tullio _ := a⃗[i] .* (I₁ .- 3) .^ i)
    ∂W∂I₂(I₂) = exp(@tullio _ := b⃗[i] .* log(I₂) .^ i)
    quadgk(∂W∂I₁, 3, I₁(λ⃗))[1] + quadgk(∂W∂I₂, 3, I₂(λ⃗))[1]
end

function ContinuumModels.StrainEnergyDensity(ψ::LambertDianiRey, I⃗, (; a⃗, b⃗), I::InvariantForm)
    ∂W∂I₁(I₁) = exp(@tullio _ := a⃗[i] .* (I₁ .- 3) .^ i)
    ∂W∂I₂(I₂) = exp(@tullio _ := b⃗[i] .* log(I₂) .^ i)
    quadgk(∂W∂I₁, 3, I⃗[1])[1] + quadgk(∂W∂I₂, 3, I⃗[2])[1]
end


function NominalStressFunction(ψ::LambertDianiRey, λ⃗::AbstractVector, (; a⃗, b⃗))
    ∂W∂I₁ = exp(@tullio _ := a⃗[i] .* (I₁(λ⃗) .- 3) .^ i)
    ∂W∂I₂ = exp(@tullio _ := b⃗[i] .* log(I₂(λ⃗)) .^ i)
    𝐒 = 2 * (I * ∂W∂I₁ - diagm(λ⃗ .^ 2)^(-2) * ∂W∂I₂)
    sᵢ = diag(𝐒)
    sᵢ = sᵢ .- sᵢ[3] .* λ⃗[3] ./ λ⃗
    return sᵢ
end

function true_stress(ψ::LambertDianiRey, λ⃗::AbstractVector, (; a⃗, b⃗))
    s(λ⃗) = NominalStressFunction(ψ, λ⃗, (a⃗=a⃗, b⃗=b⃗))
    σᵢ = map(λ⃗ᵢ -> λ⃗ᵢ .* s(λ⃗ᵢ), λ⃗)
    return σᵢ
end

function parameters(ψ::LambertDianiRey)
    return (:a⃗, :b⃗)
end

"""
Hoss Marczak I [^1]

Parameters: α, β, μ, b, n

Model: ``\\frac{\\alpha}{\\beta}(1-\\exp{-\\beta(I_1-3)})+\\frac{\\mu}{2b}\\bigg((1+\\frac{b}{n}(I_1-3))^n -1\\bigg)``

* Note: The authors suggested this model for low strains.

[^1]: > Hoss L, Marczak RJ. A new constitutive model for rubber-like materials. Mecánica Computacional. 2010;29(28):2759-73.
"""
struct HossMarczakI <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::HossMarczakI, λ⃗::AbstractVector, (; α, β, μ, b, n))
    α / β * (1 - exp(-β * (I₁(λ⃗) - 3))) + μ / (2b) * ((1 + b / n * (I₁(λ⃗) - 3))^n - 1)
end

function ContinuumModels.StrainEnergyDensity(ψ::HossMarczakI, I⃗, (; α, β, μ, b, n), I::InvariantForm)
    α / β * (1 - exp(-β * (I⃗[1] - 3))) + μ / (2b) * ((1 + b / n * (I⃗[1] - 3))^n - 1)
end

function parameters(ψ::HossMarczakI)
    return (:α, :β, :μ, :b, :n)
end

function parameter_bounds(ψ::HossMarczakI, data::AbstractHyperelasticData)
    lb = (α=-Inf, β=0, μ=-Inf, b=0, n=0)
    ub = nothing
    return (lb=lb, ub=ub)
end

"""
Hoss Marczak II [^1]

Parameters: α, β, μ, b, n, C2

Model: ``\\frac{\\alpha}{\\beta}(1-\\exp{-\\beta(I_1-3)})+\\frac{\\mu}{2b}\\bigg((1+\\frac{b}{n}(I_1-3))^n -1\\bigg)+C_2\\log(\\frac{I_2}{3})``

* Note: The authors suggests this model for high strains.

[^1]: > Hoss L, Marczak RJ. A new constitutive model for rubber-like materials. Mecánica Computacional. 2010;29(28):2759-73.
"""
struct HossMarczakII <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::HossMarczakII, λ⃗::AbstractVector, (; α, β, μ, b, n, C2))
    α / β * (1 - exp(-β * (I₁(λ⃗) - 3))) + μ / (2b) * ((1 + b / n * (I₁(λ⃗) - 3))^n - 1) + C2 * log(I₂(λ⃗) / 3)
end

function ContinuumModels.StrainEnergyDensity(ψ::HossMarczakII, I⃗, (; α, β, μ, b, n, C2), I::InvariantForm)
    α / β * (1 - exp(-β * (I⃗[1] - 3))) + μ / (2b) * ((1 + b / n * (I⃗[1] - 3))^n - 1) + C2 * log(I⃗[2] / 3)
end

function parameters(ψ::HossMarczakII)
    return (:α, :β, :μ, :b, :n, :C2)
end

function parameter_bounds(ψ::HossMarczakII, data::AbstractHyperelasticData)
    lb = (α=-Inf, β=0, μ=-Inf, b=0, n=0, C2=-Inf)
    ub = nothing
    return (lb=lb, ub=ub)
end


"""
Exp-Ln [^1]

Parameters: A, a, b

Model: ``A\\bigg[\\frac{1}{a}\\exp{(a(I_1-3))}+b(I_1-2)(1-\\log{I_1-2})-\\frac{1}{a}-b\\bigg]``

[^1]: > Khajehsaeid H, Arghavani J, Naghdabadi R. A hyperelastic constitutive model for rubber-like materials. European Journal of Mechanics-A/Solids. 2013 Mar 1;38:144-51.
"""
struct ExpLn <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::ExpLn, λ⃗::AbstractVector, (; A, a, b))
    A * (1 / a * exp(a * (I₁(λ⃗) - 3)) + b * (I₁(λ⃗) - 2) * (1 - log(I₁(λ⃗) - 2)) - 1 / a - b)
end

function ContinuumModels.StrainEnergyDensity(ψ::ExpLn, I⃗, (; A, a, b), I::InvariantForm)
    A * (1 / a * exp(a * (I⃗[1] - 3)) + b * (I⃗[1] - 2) * (1 - log(I⃗[1] - 2)) - 1 / a - b)
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
struct VanDerWaals <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::VanDerWaals, λ⃗::AbstractVector, (; μ, λm, β, α))
    I = β * I₁(λ⃗) + (1 - β) * I₂(λ⃗)
    θ = (I - 3) / (λm^2 - 3)
    μ * (-(λm^2 - 3) * log(1 - θ) + θ) - 2 / 3 * α * ((I - 3) / 2)^(3 / 2)
end

function ContinuumModels.StrainEnergyDensity(ψ::VanDerWaals, I⃗, (; μ, λm, β, α), I::InvariantForm)
    I = β * I⃗[1] + (1 - β) * I⃗[2]
    θ = (I - 3) / (λm^2 - 3)
    μ * (-(λm^2 - 3) * log(1 - θ) + θ) - 2 / 3 * α * ((I - 3) / 2)^(3 / 2)
end

function parameter_bounds(ψ::VanDerWaals, data::AbstractHyperelasticData)
    lb = (μ = 0.0, λm = sqrt(3), β = 0.0, α = 0.0)
    ub = (μ = Inf, λm = Inf, β = 1.0, α = Inf)
    return (ub = ub, lb = lb)
end

function parameters(ψ::VanDerWaals)
    return (:μ, :λm, :β, :α)
end

function constraints(ψ::VanDerWaals, data::AbstractHyperelasticData)
    I₁_max = maximum(I₁.(data.λ⃗))
    I₂_max = maximum(I₂.(data.λ⃗))
    return f(u, p) = [1 - (u.β * I₁_max + (1 - u.β) * I₂_max - 3) / (u.λm^2 - 3)]
end

"""
Gent [^1]

Parameters: μ, Jₘ

Model: ``-\\frac{\\mu J_m}{2}\\log{\\bigg(1-\\frac{I_1-3}{J_m}\\bigg)}``

[^1]: > Gent AN. A new constitutive relation for rubber. Rubber chemistry and technology. 1996 Mar;69(1):59-61.
"""
struct Gent <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Gent, λ⃗::AbstractVector, (; μ, Jₘ))
    -(μ * Jₘ) / 2 * log(1 - (I₁(λ⃗) - 3) / Jₘ)
end

function ContinuumModels.StrainEnergyDensity(ψ::Gent, I⃗, (; μ, Jₘ), I::InvariantForm)
    -(μ * Jₘ) / 2 * log(1 - (I⃗[1] - 3) / Jₘ)
end

function parameters(ψ::Gent)
    return (:μ, :Jₘ)
end

function parameter_bounds(ψ::Gent, data::AbstractHyperelasticData)
    I₁_max = maximum(I₁.(collect.(data.λ⃗)))
    Jₘ_min = I₁_max - 3
    lb = (μ=0, Jₘ=Jₘ_min)
    ub = nothing
    return (lb=lb, ub=ub)
end

"""
Takamizawa-Hayashi [^1]
From: A description of arterial wall mechanics using limiting chain extensibility constitutitive models by Horgan and Saccomandi

Parameters: c, Jₘ

Model: ``-c\\log{1-\\big(\\frac{I_1-3}{J_m}\\big)^2}``

[^1]: > Takamizawa K, Hayashi K. Strain energy density function and uniform strain hypothesis for arterial mechanics. Journal of biomechanics. 1987 Jan 1;20(1):7-17.
"""
struct TakamizawaHayashi <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::TakamizawaHayashi, λ⃗::AbstractVector, (; c, Jₘ))
    -c * log(1 - ((I₁(λ⃗) - 3) / Jₘ)^2)
end

function ContinuumModels.StrainEnergyDensity(ψ::TakamizawaHayashi, I⃗, (; c, Jₘ), I::InvariantForm)
    -c * log(1 - ((I⃗[1] - 3) / Jₘ)^2)
end

function parameters(ψ::TakamizawaHayashi)
    return (:c, :Jₘ)
end

function parameter_bounds(ψ::TakamizawaHayashi, data::AbstractHyperelasticData)
    I₁_max = maximum(I₁.(collect.(data.λ⃗)))
    Jₘ_min = I₁_max - 3
    lb = (c=-Inf, Jₘ=Jₘ_min)
    ub = nothing
    return (lb=lb, ub=ub)
end

"""
Yeoh-Fleming [^1]

Parameters: A, B, C10, Im

Model: ``\\frac{A}{B}(1-\\exp{-B(I_1-3)}) - C_{10}(I_m-3)\\log{1-\\frac{I_1-3}{I_m-3}}``

[^1]: >  Yeoh OH, Fleming PD. A new attempt to reconcile the statistical and phenomenological theories of rubber elasticity. Journal of Polymer Science Part B: Polymer Physics. 1997 Sep 15;35(12):1919-31.
"""
struct YeohFleming <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::YeohFleming, λ⃗::AbstractVector, (; A, B, C10, Im))
    A / B * (1 - exp(-B * (I₁(λ⃗) - 3))) - C10 * (Im - 3) * log(1 - ((I₁(λ⃗) - 3) / (Im - 3)))
end

function ContinuumModels.StrainEnergyDensity(ψ::YeohFleming, I⃗, (; A, B, C10, Im), I::InvariantForm)
    A / B * (1 - exp(-B * (I⃗[1] - 3))) - C10 * (Im - 3) * log(1 - ((I⃗[1] - 3) / (Im - 3)))
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
struct PucciSaccomandi <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::PucciSaccomandi, λ⃗::AbstractVector, (; K, μ, Jₘ))
    K * log(I₂(λ⃗) / 3) - μ * Jₘ / 2 * log(1 - (I₁(λ⃗) - 3) / Jₘ)
end

function ContinuumModels.StrainEnergyDensity(ψ::PucciSaccomandi, I⃗, (; K, μ, Jₘ), I::InvariantForm)
    K * log(I⃗[2] / 3) - μ * Jₘ / 2 * log(1 - (I⃗[1] - 3) / Jₘ)
end

function parameters(ψ::PucciSaccomandi)
    return (:K, :μ, :Jₘ)
end

function parameter_bounds(ψ::PucciSaccomandi, data::AbstractHyperelasticData)
    I₁_max = maximum(I₁.(collect.(data.λ⃗)))
    Jₘ_min = I₁_max - 3
    lb = (K=-Inf, μ=-Inf, Jₘ=Jₘ_min)
    ub = nothing
    return (lb=lb, ub=ub)
end

"""
Horgan Saccomandi Model [^1][^2]

Parameters: μ, J

Model: ``-\\frac{\\mu J}{2}\\log\\bigg(\\frac{J^3-J^2I_1+JI_2-1}{(J-1)^3}\\bigg)``

[^1]: > Horgan CO, Saccomandi G. Constitutive models for compressible nonlinearly elastic materials with limiting chain extensibility. Journal of Elasticity. 2004 Nov;77(2):123-38.\
[^2]: > Horgan CO, Saccomandi G. Constitutive models for atactic elastomers. InWaves And Stability In Continuous Media 2004 (pp. 281-294).
"""
struct HorganSaccomandi <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::HorganSaccomandi, λ⃗::AbstractVector, (; μ, J))
    -μ * J / 2 * log((J^3 - J^2 * I₁(λ⃗) + J * I₂(λ⃗) - 1) / (J - 1)^3)
end

function ContinuumModels.StrainEnergyDensity(ψ::HorganSaccomandi, I⃗, (; μ, J), I::InvariantForm)
    -μ * J / 2 * log((J^3 - J^2 * I⃗[1] + J * I⃗[2] - 1) / (J - 1)^3)
end

function parameters(ψ::HorganSaccomandi)
    return (:μ, :J)
end

function parameter_bounds(ψ::HorganSaccomandi, data::AbstractHyperelasticData)
    J_min = maximum(maximum.(map(x->x.^2,data.λ⃗)))
    lb = (μ=-Inf, J=J_min)
    ub = nothing
    return (lb=lb, ub=ub)
end

# function constraints(ψ::HorganSaccomandi, data::AbstractHyperelasticData)
#     I₁_max = maximum(I₁.(data.λ⃗))
#     I₂_max = maximum(I₂.(data.λ⃗))
#     f(u, p) = [(u.J^3 - u.J^2 * I₁_max + u.J * I₂_max - 1) / (u.J - 1)^3]
#     return f
# end


"""
Beatty Model [^1]

Parameters: G₀, Iₘ

Model: ``-\\frac{G_0 I_m(I_m-3)}{2(2I_m-3)}\\log\\bigg(\\frac{1-\\frac{I_1-3}{I_m-3}}{1+\\frac{I_1-3}{I_m}} \\bigg)``

[^1]: > Beatty MF. On constitutive models for limited elastic, molecular based materials. Mathematics and mechanics of solids. 2008 Jul;13(5):375-87.
"""
struct Beatty <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Beatty, λ⃗::AbstractVector, (; G₀, Iₘ))
    -G₀ * Iₘ * (Iₘ - 3) / 2 / (2Iₘ - 3) * log((1 - (I₁(λ⃗) - 3) / (Iₘ - 3)) / (1 + (I₁(λ⃗) - 3) / (Iₘ)))
end

function ContinuumModels.StrainEnergyDensity(ψ::Beatty, I⃗, (; G₀, Iₘ), I::InvariantForm)
    -G₀ * Iₘ * (Iₘ - 3) / 2 / (2Iₘ - 3) * log((1 - (I⃗[1] - 3) / (Iₘ - 3)) / (1 + (I⃗[1] - 3) / (Iₘ)))
end

function parameters(ψ::Beatty)
    return (:G₀, :Iₘ)
end

"""
Horgan Murphy Model [^1]

Parameters: μ, Jₘ, c

Model: ``-\\frac{2\\mu J_m}{c^2}\\log\\bigg(1-\\frac{\\lambda_1^c+\\lambda_2^c+\\lambda_3^c-3}{J_m})``

[^1]: > Horgan CO, Murphy JG. Limiting chain extensibility constitutive models of Valanis–Landel type. Journal of Elasticity. 2007 Feb;86(2):101-11.
"""
struct HorganMurphy <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::HorganMurphy, λ⃗::AbstractVector, (; μ, Jₘ, c))
    -2 * μ * Jₘ / c^2 * log(1 - (sum(λ⃗ .^ c) - 3) / Jₘ)
end

function parameters(ψ::HorganMurphy)
    return (:μ, :Jₘ, :c)
end

function parameter_bounds(ψ::HorganMurphy, data::AbstractHyperelasticData)
    lb = (μ=-Inf, Jₘ=0, c=0)
    ub = (μ=Inf, Jₘ=Inf, c=Inf)
    return (lb=lb, ub=ub)
end

function constraints(ψ::HorganMurphy, data::AbstractHyperelasticData)
    function f(u, p)
        max_sum = maximum(λ⃗ -> sum(λ⃗ .^ u.c), data.λ⃗)
        [1 - (max_sum - 3) / u.Jₘ]
    end
    return f
end

"""
Valanis-Landel [^1]

Parameters: μ

Model: ``2\\mu\\sum\\limits_{1}^{3}(\\lambda_i(\\log\\lambda_i -1))``

[^1]: Valanis KC, Landel RF. The strain‐energy function of a hyperelastic material in terms of the extension ratios. Journal of Applied Physics. 1967 Jun;38(7):2997-3002.
"""
struct ValanisLandel <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::ValanisLandel, λ⃗::AbstractVector, (; μ))
    2 * μ * sum(λ⃗ .* (log.(λ⃗) .- 1))
end

function parameters(ψ::ValanisLandel)
    return (:μ,)
end

"""
Peng - Landel [^1]

Parameters: E

Model: ``E\\sum\\limits_{i=1}^{3}\\bigg[\\lambda_i - 1 - \\log(\\lambda_i) - \\frac{1}{6}\\log(\\lambda_i)^2 + \\frac{1}{18}\\log(\\lambda_i)^3-\\frac{1}{216}\\log(\\lambda_i)^4\\bigg]``

[^1]: > Peng TJ, Landel RF. Stored energy function of rubberlike materials derived from simple tensile data. Journal of Applied Physics. 1972 Jul;43(7):3064-7.
"""
struct PengLandel <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::PengLandel, λ⃗::AbstractVector, (; E))
    @tullio _ := (λ⃗[i] - 1 - log(λ⃗[i]) - 1 / 6 * log(λ⃗[i])^2 + 1 / 18 * log(λ⃗[i])^3 - 1 / 216 * log(λ⃗[i])^4) * E
end

function parameters(ψ::PengLandel)
    return (:E,)
end

"""
Ogden [^1]

Parameters: μ⃗, α⃗

Model: ``\\sum\\limits_{i=1}^{N}\\frac{\\mu_i}{\\alpha_i}(\\lambda_1^{\\alpha_i}+\\lambda_2^{\\alpha_i}+\\lambda_3^{\\alpha_i}-3)``

[^1]: > Ogden RW. Large deformation isotropic elasticity–on the correlation of theory and experiment for incompressible rubberlike solids. Proceedings of the Royal Society of London. A. Mathematical and Physical Sciences. 1972 Feb 1;326(1567):565-84.
"""
struct Ogden <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Ogden, λ⃗::AbstractVector, (; μ⃗, α⃗))
    @tullio _ := μ⃗[i] / α⃗[i] * (sum(λ⃗ .^ α⃗[i]) - 3)
end

function parameters(ψ::Ogden)
    return (:μ⃗, :α⃗)
end

"""
Attard [^1]

Parameters: A⃗, B⃗

Model: ``\\sum\\limits_{i=1}^N\\frac{A_i}{2i}(\\lambda_1^{2i}+\\lambda_2^{2i}+\\lambda_3^{2i}-3) + \\frac{B_i}{2i}(\\lambda_1^{-2i}+\\lambda_2^{-2i}+\\lambda_3^{-2i}-3)``

[^1]: > Attard MM, Hunt GW. Hyperelastic constitutive modeling under finite strain. International Journal of Solids and Structures. 2004 Sep 1;41(18-19):5327-50.
"""
struct Attard <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Attard, λ⃗::AbstractVector, (; A⃗, B⃗))
    @assert length(A⃗) == length(B⃗) "Length of A and B are not equal"
    @tullio _ := A⃗[i] / 2 / i * (sum(λ⃗ .^ (2i)) - 3) + B⃗[i] / 2 / i * (sum(λ⃗ .^ (-2i)) - 3)
end

function parameters(ψ::Attard)
    return (:A⃗, :B⃗)
end

"""
Shariff [^1]

Parameters: E, α⃗

Model:
``E\\sum\\limits_{i=1}^3\\sum\\limits_{j=1}^{N}\\alpha_j \\Phi_j(\\lambda_i)``

[^1]: > Shariff MH. Strain energy function for filled and unfilled rubberlike material. Rubber chemistry and technology. 2000 Mar;73(1):1-8.
"""
struct Shariff <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Shariff, λ⃗::AbstractVector, (; E, α⃗))
    ϕ = []
    c(j, r) = factorial(j) / factorial(r) / factorial(j - r)
    for j in eachindex(α⃗)
        if j == 0
            push!(ϕ, x -> log(x)^2 / 3)
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
    E * (@tullio _ := ϕ[i](λ⃗[j]).*α⃗[i])
end

function parameters(ψ::Shariff)
    return (:E, :α⃗)
end

"""
Arman - Narooei [^1]

Parameters: A⃗, B⃗, m⃗, n⃗, α⃗, β⃗

Model: ``\\sum\\limits_{i=1}^{N} A_i\\big[\\exp{m_i(\\lambda_1^{\\alpha_i}+\\lambda_2^{\\alpha_i}+\\lambda_3^{\\alpha_i}-3)}-1] + B_i\\big[\\exp{n_i(\\lambda_1^{-\\beta_i}+\\lambda_2^{-\\beta_i}+\\lambda_3^{-\\beta_i}-3)}-1]``

[^1]: > Narooei K, Arman M. Modification of exponential based hyperelastic strain energy to consider free stress initial configuration and Constitutive modeling. Journal of Computational Applied Mechanics. 2018 Jun 1;49(1):189-96.
"""
struct ArmanNarooei <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::ArmanNarooei, λ⃗::AbstractVector, (; A⃗, B⃗, m⃗, n⃗, α⃗, β⃗))
    @assert length(A⃗) == length(B⃗) == length(m⃗) == length(n⃗) == length(α⃗) == length(β⃗) "Length of A, B, m, n, α and β are not equal"
    @tullio _ := A⃗[i] * (exp(m⃗[i] * (sum(λ⃗ .^ α⃗[i]) - 3)) - 1) + B⃗[i] * (exp(n⃗[i] * (sum(λ⃗ .^ (-β⃗[i])) - 3)) - 1)
end

function parameters(ψ::ArmanNarooei)
    return (:A⃗, :B⃗, :m⃗, :n⃗, :α⃗, :β⃗)
end

"""
Continuum Hybrid [^1]

Parameters: K₁, K₂, α, μ

Model: ``K_1(I_1-3)+K_2\\log\\frac{I_2}{3}+\\frac{\\mu}{\\alpha}(\\lambda_1^\\alpha+\\lambda_2^\\alpha+\\lambda^\\alpha-3)``

[^1]: > Beda T, Chevalier Y. Hybrid continuum model for large elastic deformation of rubber. Journal of applied physics. 2003 Aug 15;94(4):2701-6.
"""
struct ContinuumHybrid <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::ContinuumHybrid, λ⃗::AbstractVector, (; K₁, K₂, α, μ))
    K₁ * (I₁(λ⃗) - 3) + K₂ * log(I₂(λ⃗) / 3) + μ / α * (sum(λ⃗ .^ α) - 3)
end

function parameters(ψ::ContinuumHybrid)
    return (:K₁, :K₂, :α, :μ)
end

"""
Bechir-4 Term [^1]

Parameters: C11, C12, C21, C22

Model: ``C_1^1(I_1-3)+\\sum\\limits_{n=1}^{2}\\sum\\limits_{r=1}^{2}C_n^{r}(\\lambda_1^{2n}+\\lambda_2^{2n}+\\lambda_3^{2n}-3)^r``

[^1]: > Khajehsaeid H, Arghavani J, Naghdabadi R. A hyperelastic constitutive model for rubber-like materials. European Journal of Mechanics-A/Solids. 2013 Mar 1;38:144-51.
"""
struct Bechir4Term <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Bechir4Term, λ⃗::AbstractVector, (; C11, C12, C21, C22))
    C = [C11 C12; C21 C22]
    C[1, 1] * (I₁(λ⃗) - 3) + sum(n -> sum(r -> C[n, r] * (sum(λ⃗ .^ (2n))), 1:2), 1:2)
end

function parameters(ψ::Bechir4Term)
    return (:C11, :C12, :C21, :C22)
end

"""
Constrained Junction [^1][^2]

Parameters: Gc, νkT, κ

Model: ``G_c (I_1-3)+ \\frac{\\nu k T}{2}(\\sum\\limits_{i=1}^{3}\\kappa\\frac{\\lambda_i-1}{\\lambda_i^2+\\kappa}+\\log{\\frac{\\lambda_i^2+\\kappa}{1+\\kappa}}-\\log{\\lambda_i^2})``

[^1]: > Flory PJ, Erman B. Theory of elasticity of polymer networks. 3. Macromolecules. 1982 May;15(3):800-6.
[^2]: > Erman B, Flory PJ. Relationships between stress, strain, and molecular constitution of polymer networks. Comparison of theory with experiments. Macromolecules. 1982 May;15(3):806-11.
"""
struct ConstrainedJunction <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::ConstrainedJunction, λ⃗::AbstractVector, (; Gc, μkT, κ))
    Gc * (I₁(λ⃗) - 3) + μkT / 2 * sum(i -> κ * (λ⃗[i] - 1) / (λ⃗[i]^2 + κ) + log((λ⃗[i]^2 + κ) / (1 + κ)) - log(λ⃗[i]^2), 1:3)
end

function parameters(ψ::ConstrainedJunction)
    return (:Gc, :μkT, :κ)
end

function parameter_bounds(ψ::ConstrainedJunction, data::AbstractHyperelasticData)
    λ_min = minimum(minimum.(collect.(data.λ⃗)))
    κ_min = -λ_min^2
    lb = (Gc=-Inf, μkT=-Inf, κ=κ_min)
    ub = nothing
    return (lb=lb, ub=ub)
end
"""
Edward-Vilgis [^1]

Parameters: Ns, Nc, α, η

Model: ``\\frac{1}{2}N_C\\Bigg[\\frac{(1-\\alpha^2)I_1}{1-\\alpha^2I_1}+\\log(1-\\alpha^2I_1)\\Bigg]+\\frac{1}{2}N_S\\Bigg[\\sum_{i=1}^{3}\\Big\\{\\frac{(1+\\eta)(1-\\alpha^2)\\lambda_i^2}{( 1+\\eta\\lambda_i^2)(1-\\alpha^2I_1)}+\\log(1+\\eta\\lambda_i^2)\\Big\\}+\\log(1-\\alpha^2I_1)\\Bigg]``

[^1]: > Edwards SF, Vilgis T. The effect of entanglements in rubber elasticity. Polymer. 1986 Apr 1;27(4):483-92.
"""
struct EdwardVilgis <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::EdwardVilgis, λ⃗::AbstractVector, (; Ns, Nc, α, η))
    0.5 * Nc * ((1 - α^2) * I₁(λ⃗) / (1 - α^2 * I₁(λ⃗)) + log(1 - α^2 * I₁(λ⃗))) + 0.5 * Ns * ((1 + η) * (1 - α^2) * λ⃗[1] / (1 + η * λ⃗[1]^2) / (1 - α^2 * I₁(λ⃗)) + log(1 + η * λ⃗[1]^2) + (1 + η) * (1 - α^2) * λ⃗[2] / (1 + η * λ⃗[2]^2) / (1 - α^2 * I₁(λ⃗)) + log(1 + η * λ⃗[2]^2) + (1 + η) * (1 - α^2) * λ⃗[3] / (1 + η * λ⃗[3]^2) / (1 - α^2 * I₁(λ⃗)) + log(1 + η * λ⃗[3]^2) + log(1 - α^2 * I₁(λ⃗)))
end

function parameters(ψ::EdwardVilgis)
    return (:Ns, :Nc, :α, :η)
end

function parameter_bounds(ψ::EdwardVilgis, data::AbstractHyperelasticData)
    I₁_max = maximum(I₁.(collect.(data.λ⃗)))
    λ_max = maximum(maximum.(data.λ⃗))
    η_min = -1 / λ_max^2
    α_max = sqrt(1 / I₁_max)
    lb = (Ns=-Inf, Nc=-Inf, α=-α_max, η=η_min)
    ub = (Ns=Inf, Nc=Inf, α=α_max, η=Inf)
    return (lb=lb, ub=ub)
end
"""
MCC (modified constrained chain) [^1]

Parameters:

Model:

``\\frac{1}{2}\\zeta k T \\sum\\limits_{i=1}^{3}(\\lambda_i^2-1)+\\frac{1}{2}\\mu k T\\sum\\limits_{i=1}^{3}[B_i+D_i-\\log{(1+B_i)}-\\log{(1+D_i)}]``

``B_i = \\frac{\\kappa^2(\\lambda_i^2-1)}{(\\lambda_i^2+\\kappa)^2}``

``D_i = \\frac{\\lambda_i^2 B_i}{\\kappa}``

[^1]: > Erman B, Monnerie L. Theory of elasticity of amorphous networks: effect of constraints along chains. Macromolecules. 1989 Aug;22(8):3342-8.
"""
struct MCC <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::MCC, λ⃗::AbstractVector, (; ζkT, μkT, κ))
    @tullio B[i] := κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2)
    @tullio D[i] := λ⃗[i]^2 * B[i] / κ
    @tullio W1 := λ⃗[i]^2 - 1
    @tullio W2 := B[i] - log(1 + B[i])
    @tullio W3 := D[i] - log(1 + D[i])
    return 1 / 2 * ζkT * W1 + 1 / 2 * μkT * (W2 + W3)
    # W(λ⃗) = 1 / 2 * ζkT * sum(i -> λ⃗[i]^2 - 1, 1:3) + 1 / 2 * μkT * sum(i -> κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2) + (λ⃗[i]^2 * (κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2)) / κ) - log(1 + (κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2))) - log(1 + (λ⃗[i]^2 * (κ^2 * (λ⃗[i]^2 - 1) * (λ⃗[i]^2 + κ)^(-2)) / κ)), 1:3)
end

function parameters(ψ::MCC)
    return (:ζkT, :μkT, :κ)
end

function parameter_bounds(ψ::MCC, data::AbstractHyperelasticData)
    lb = (ζkT=-Inf, μkT=-Inf, κ=0)
    ub = nothing
    return (lb=lb, ub=ub)
end

"""
Tube [^1]

Parameters: Gc, Ge, β

Model: ``\\sum\\limits_{i=1}^{3}\\frac{G_c}{2}(\\lambda_i^2-1)+\\frac{2Ge}{\\beta^2}(\\lambda_i^{-\\beta}-1)``

[^1]: > Heinrich G, Kaliske M. Theoretical and numerical formulation of a molecular based constitutive tube-model of rubber elasticity. Computational and Theoretical Polymer Science. 1997 Jan 1;7(3-4):227-41.
"""
struct Tube <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::Tube, λ⃗::AbstractVector, (; Gc, Ge, β))
    @tullio _ := Gc / 2 * (λ⃗[i]^2 - 1) + 2Ge / β^2 * (λ⃗[i]^(-β) - 1)
end

function parameters(ψ::Tube)
    return (:Gc, :Ge, :β)
end

"""
Nonaffine - Tube [^1]

Parameters: Gc, Ge

Model: ``G_c \\sum\\limits_{i=1}^{3}\\frac{\\lambda_i^2}{2}+G_e\\sum\\limits_{i=1}^{3}\\lambda_i+\\frac{1}{\\lambda_i}``

[^1]: > Rubinstein M, Panyukov S. Nonaffine deformation and elasticity of polymer networks. Macromolecules. 1997 Dec 15;30(25):8036-44.
"""
struct NonaffineTube <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::NonaffineTube, λ⃗::AbstractVector, (; Gc, Ge))
    Gc * sum(λ⃗ .^ 2 ./ 2) + Ge * sum(λ⃗ .+ 1 ./ λ⃗)
end

function parameters(ψ::NonaffineTube)
    return (:Gc, :Ge)
end

"""
Three Chain Model [^1]

Note: The field `ℒinv` can be set to change the inverse Langevin function approximation used. Currently, the default choice is the Pade 3/2 Approximation from Cohen 1991 [^2]

Parameters: μ, N

Model: `` \\frac{\\mu\\sqrt{N}}{3}\\sum\\limits_{i=1}^{3}\\bigg(\\lambda_i\\beta_i+\\sqrt{N}\\log\\bigg(\\frac{\\beta_i}{\\sinh \\beta_i}\\bigg)\\bigg)``

[^1]: > James HM, Guth E. Theory of the elastic properties of rubber. The Journal of Chemical Physics. 1943 Oct;11(10):455-81.
[^2]: > Cohen A. A Padé approximant to the inverse Langevin function. Rheologica acta. 1991 May;30(3):270-3.

"""
struct ThreeChainModel <: AbstractHyperelasticModel
    ℒinv::Function
    ThreeChainModel(; ℒinv::Function=TreloarApproximation) = new(ℒinv)
end

function ContinuumModels.StrainEnergyDensity(ψ::ThreeChainModel, λ⃗::AbstractVector, (; μ, N))
    μ * sqrt(N) / 3 * sum(λ⃗ .* ψ.ℒinv.(λ⃗ ./ sqrt(N)) .+ sqrt(N) .* log.((ψ.ℒinv.(λ⃗ ./ sqrt(N))) ./ (sinh.(ψ.ℒinv.(λ⃗ ./ sqrt(N))))))
end

function parameters(ψ::ThreeChainModel)
    return (:μ, :N)
end

function parameter_bounds(ψ::ThreeChainModel, data::AbstractHyperelasticData)
    λ_max = maximum(maximum.(collect.(data.λ⃗)))
    N_min = λ_max^2
    lb = (μ=-Inf, N=N_min)
    ub = nothing
    return (lb=lb, ub=ub)
end

"""
Arruda Boyce [^1]

Note: The field `ℒinv` can be set to change the inverse Langevin function approximation used. Currently, the default choice is the Pade 3/2 Approximation from Cohen 1991 [^2]

Parameters: μ, N

Model: ``\\mu\\bigg(\\frac{1}{2}(I_1-3)+\\frac{I_1^2-9}{20N}+\\frac{11(I_1^3-27)}{1050N^2}+\\frac{19(I_1^4-81)}{7000N^3}+\\frac{519(I_1^5-243)}{673750N^4}\\bigg)``


[^1]: > Arruda EM, Boyce MC. A three-dimensional constitutive model for the large stretch behavior of rubber elastic materials. Journal of the Mechanics and Physics of Solids. 1993 Feb 1;41(2):389-412.

[^2]: > Cohen A. A Padé approximant to the inverse Langevin function. Rheologica acta. 1991 May;30(3):270-3.
"""
struct ArrudaBoyce <: AbstractHyperelasticModel
    ℒinv::Function
    ArrudaBoyce(; ℒinv::Function=TreloarApproximation) = new(ℒinv)
end

function ContinuumModels.StrainEnergyDensity(ψ::ArrudaBoyce, λ⃗::AbstractVector, (; μ, N))
    rchain_Nl = √(I₁(λ⃗) / 3 / N)
    β = ψ.ℒinv(rchain_Nl)
    μ * N * (rchain_Nl * β + log(β / sinh(β)))
end

# function true_stress(ψ::ArrudaBoyce, λ⃗, (; μ, N))
#     λch = sqrt(I₁(λ⃗)/3)
#     @tullio σ[i] := μ * √(N) * (λ⃗[i]^2 - λch^2)/(λch)*ψ.ℒinv(λch/sqrt(N))
#     return σ
# end

# function ContinuumModels.StrainEnergyDensity(ψ::ArrudaBoyce, I⃗, (; μ, N), I::InvariantForm)
#     rchain_Nl = √(I⃗[1] / 3 / N)
#     β = ψ.ℒinv(rchain_Nl)
#     μ * N * (rchain_Nl * β + log(β / sinh(β)))
# end

function parameters(ψ::ArrudaBoyce)
    return (:μ, :N)
end

function parameter_bounds(ψ::ArrudaBoyce, data::AbstractHyperelasticData)
    I₁_max = maximum(I₁.(collect.(data.λ⃗)))
    N_max = 11 / 35 * I₁_max # old
    N_max = I₁_max / 3
    lb = (μ=-Inf, N=N_max)
    ub = nothing
    return (lb=lb, ub=ub)
end

"""
Modified Flory Erman [^1]

Parameters: μ, N, κ

Model: ``W_{\\text{Arruda-Boyce}}+\\sum\\limits_{i=1}^{3}\\frac{\\mu}{2}[B_i+D_i]

[^1]: > Edwards SF. The statistical mechanics of polymerized material. Proceedings of the Physical Society (1958-1967). 1967 Sep 1;92(1):9.
"""
struct ModifiedFloryErman <: AbstractHyperelasticModel
    ℒinv::Function
    ModifiedFloryErman(; ℒinv::Function=TreloarApproximation) = new(ℒinv)
end

function ContinuumModels.StrainEnergyDensity(ψ::ModifiedFloryErman, λ⃗::AbstractVector, (; μ, N, κ))
    WAB = ContinuumModels.StrainEnergyDensity(ArrudaBoyce(ℒinv=ψ.ℒinv), λ⃗, (μ=μ, N=N))
    @tullio B[i] := κ^2 * (λ⃗[i]^2 - 1) / (λ⃗[i]^2 + κ)^2
    @tullio D[i] := λ⃗[i]^2 * B[i] / κ
    @tullio W2 := B[i] + D[i] - log(B[i] + 1) - log(D[i] + 1)
    WAB + W2
end

function parameters(ψ::ModifiedFloryErman)
    return (:μ, :N, :κ)
end

function parameter_bounds(ψ::ModifiedFloryErman, data::AbstractHyperelasticData)
    I₁_max = maximum(I₁.(collect.(data.λ⃗)))
    N_max = 11 / 35 * I₁_max # old
    N_max = I₁_max / 3
    lb = (μ=-Inf, N=N_max, κ=0)
    ub = nothing
    return (lb=lb, ub=ub)
end

"""
Extended Tube Model [^1]

Parameters: Gc, Ge, δ, β

Model: ``\\frac{G_c}{2}\\bigg[\\frac{(1-\\delta^2)(I_1-3)}{1-\\delta^2(I_1-3)}+\\log{(1-\\delta^2(I_1-3))}\\bigg]+\\frac{2G_e}{\\beta^2}\\sum\\limits_{i=1}^{3}(\\lambda_i^{-\\beta}-1)``

[^1]: > Kaliske M, Heinrich G. An extended tube-model for rubber elasticity: statistical-mechanical theory and finite element implementation. Rubber Chemistry and Technology. 1999 Sep;72(4):602-32.
"""
struct ExtendedTubeModel <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::ExtendedTubeModel, λ⃗::AbstractVector, (; Gc, Ge, δ, β))
    Gc / 2 * ((1 - δ^2) * (I₁(λ⃗) - 3) / (1 - δ^2 * (I₁(λ⃗) - 3)) + log(1 - δ^2 * (I₁(λ⃗) - 3))) + 2 * Ge / β^2 * sum(λ⃗ .^ (-β) .- 1)
end

function parameters(ψ::ExtendedTubeModel)
    return (:Gc, :Ge, :δ, :β)
end

function parameter_bounds(ψ::ExtendedTubeModel, data::AbstractHyperelasticData)
    I₁_max = maximum(I₁.(collect.(data.λ⃗)))
    δ_max = sqrt(1 / (I₁_max - 3))
    lb = (Gc=-Inf, Ge=-Inf, δ=-δ_max, β=0)
    ub = (Gc=Inf, Ge=Inf, δ=δ_max, β=Inf)
    return (lb=lb, ub=ub)
end

"""
ABGI [^1][^2]

Parameters: μ, N, Ge, n

Model: ``W_{Arruda-Boyce} + G_e\\frac{\\lambda_1^n+\\lambda_2^2+\\lambda_3^2-3}{n}``

[^1]: > Meissner B, Matějka L. A Langevin-elasticity-theory-based constitutive equation for rubberlike networks and its comparison with biaxial stress–strain data. Part I. Polymer. 2003 Jul 1;44(16):4599-610.
[^2]: > Meissner B, Matějka L. A Langevin-elasticity-theory-based constitutive equation for rubberlike networks and its comparison with biaxial stress–strain data. Part I. Polymer. 2003 Jul 1;44(16):4599-610.
"""
struct ABGI <: AbstractHyperelasticModel
    ℒinv::Function
    ABGI(; ℒinv::Function=TreloarApproximation) = new(ℒinv)
end

function ContinuumModels.StrainEnergyDensity(ψ::ABGI, λ⃗::AbstractVector, (; μ, N, Ge, n))
    WAB = ContinuumModels.StrainEnergyDensity(ArrudaBoyce(ℒinv=ψ.ℒinv), λ⃗::AbstractVector, (μ=μ, N=N))
    WAB + Ge * (sum(λ⃗ .^ n) - 3) / n
end

function parameters(ψ::ABGI)
    return (:μ, :N, :Ge, :n)
end

function parameter_bounds(ψ::ABGI, data::AbstractHyperelasticData)
    I₁_max = maximum(I₁.(collect.(data.λ⃗)))
    lb = (μ=-Inf, N=11 / 35 * I₁_max, Ge=-Inf, n=0)
    ub = nothing
    return (lb=lb, ub=ub)
end
"""
Non-Affine Micro-Sphere [^1]

Note: The field `ℒinv` can be set to change the inverse Langevin function approximation used. Currently, the default choice is the Pade 3/2 Approximation from Cohen 1991 [^2]

Parameters: μ, N, p, U, q

Model: See Paper

---
[^1]: > Miehe C, Göktepe S, Lulei F. A micro-macro approach to rubber-like materials—part I: the non-affine micro-sphere model of rubber elasticity. Journal of the Mechanics and Physics of Solids. 2004 Nov 1;52(11):2617-60.
[^2]: > Cohen A. A Padé approximant to the inverse Langevin function. Rheologica acta. 1991 May;30(3):270-3.
"""
struct NonaffineMicroSphere <: AbstractHyperelasticModel
    ℒinv::Function
    NonaffineMicroSphere(; ℒinv::Function=TreloarApproximation) = new(ℒinv)
end

function ContinuumModels.StrainEnergyDensity(ψ::NonaffineMicroSphere, λ⃗::AbstractVector, (; μ, N, p, U, q))
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
    F = diagm(λ⃗)
    @tullio t⃗[i] := F * r⃗[i]
    @tullio n⃗[i] := inv(F') * r⃗[i]
    @tullio λ̄[i] := norm(t⃗[i])
    @tullio ν̄[i] := norm(n⃗[i])
    @tullio λ := (λ̄[i]^p) * w[i]# |> Base.Fix2(^, (1 / p))
    λr = λ^(1 / p) / √N
    β = ψ.ℒinv(λr)
    @tullio ν := ν̄[i]^q * w[i]# |> Base.Fix2(^, 1 / q)
    return N * U * μ * ν^(1 / q) + N * μ * (λr * β + log(β / sinh(β)))
end

function parameters(ψ::NonaffineMicroSphere)
    return (:μ, :N, :p, :U, :q)
end

function parameter_bounds(ψ::NonaffineMicroSphere, data::AbstractHyperelasticData)
    lb = (μ=-Inf, N=0, p=0, U=0, q=0)
    ub = nothing
    return (lb=lb, ub=ub)
end

"""
Affine Micro-Sphere [^1]

Note: The field `ℒinv` can be set to change the inverse Langevin function approximation used. Currently, the default choice is the Pade 3/2 Approximation from Cohen 1991 [^2]

Parameters: μ, N, p, U, q

Model: See Paper

---
[^1]: > Miehe C, Göktepe S, Lulei F. A micro-macro approach to rubber-like materials—part I: the non-affine micro-sphere model of rubber elasticity. Journal of the Mechanics and Physics of Solids. 2004 Nov 1;52(11):2617-60.
[^2]: > Cohen A. A Padé approximant to the inverse Langevin function. Rheologica acta. 1991 May;30(3):270-3.
"""
struct AffineMicroSphere <: AbstractHyperelasticModel
    ℒinv::Function
    AffineMicroSphere(; ℒinv::Function=TreloarApproximation) = new(ℒinv)
end

function ContinuumModels.StrainEnergyDensity(ψ::AffineMicroSphere, λ⃗::AbstractVector, (; μ, N, p, U, q))
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

    F = diagm(λ⃗)
    @tullio t⃗[i] := F * r⃗[i]
    @tullio n⃗[i] := inv(F') * r⃗[i]
    @tullio λ̄[i] := norm(t⃗[i])
    @tullio ν̄[i] := norm(n⃗[i])
    @tullio λ := (λ̄[i]) * w[i]# |> Base.Fix2(^, (1 / p))
    λr = λ^(1 / p) / √N
    β = ψ.ℒinv(λr)
    @tullio ν := ν̄[i]^q * w[i]# |> Base.Fix2(^, 1 / q)
    return N * U * μ * ν^(1 / q) + N * μ * (λr * β + log(β / sinh(β)))
end

function parameters(ψ::AffineMicroSphere)
    return (:μ, :N, :p, :U, :q)
end

"""
Bootstrapped 8Chain Model [^1][^2]

Note: The field `ℒinv` can be set to change the inverse Langevin function approximation used. Currently, the default choice is the Pade 3/2 Approximation from Cohen 1991 [^3]

Parameters: μ, N

Model: ``W_8(\\frac{\\sum\\lambda}{\\sqrt{3N}}-\\frac{\\lambda_{chain}}{\\sqrt{N}})+W_{8}(\\frac{\\lambda_{chain}}{\\sqrt{N}})``

``W_8(x) = \\mu N (x \\mathcal{L}^{-1}(x) + \\log\\frac{\\mathcal{L}^{-1}(x)}{\\sinh\\mathcal{L}^{-1}(x)})``

``\\lambda_{chain} = \\sqrt{\\frac{I_1}{3}}``

[^1]: > Miroshnychenko D, Green WA, Turner DM. Composite and filament models for the mechanical behaviour of elastomeric materials. Journal of the Mechanics and Physics of Solids. 2005 Apr 1;53(4):748-70.
[^2]: > Miroshnychenko D, Green WA. Heuristic search for a predictive strain-energy function in nonlinear elasticity. International Journal of Solids and Structures. 2009 Jan 15;46(2):271-86.
[^3]: >
"""
struct Bootstrapped8Chain <: AbstractHyperelasticModel
    ℒinv::Function
    Bootstrapped8Chain(; ℒinv::Function=TreloarApproximation) = new(ℒinv)
end

function ContinuumModels.StrainEnergyDensity(ψ::Bootstrapped8Chain, λ⃗::AbstractVector, (; μ, N))
    function W8(x)
        β = ψ.ℒinv(x)
        μ * N * (x * β + log(β / sinh(β)))
    end
    λchain = √(I₁(λ⃗) / 3)
    W8(sum(λ⃗) / √(3N) - λchain / √(N)) + W8(λchain / √(N))
end

function parameters(ψ::Bootstrapped8Chain)
    return (:μ, :N)
end

function parameter_bounds(ψ::Bootstrapped8Chain, data::AbstractHyperelasticData)
    I₁_max = maximum(I₁.(collect.(data.λ⃗)))
    N_min = I₁_max / 3
    lb = (μ=-Inf, N=N_min)
    ub = nothing
    return (lb=lb, ub=ub)
end

"""
Davidson - Goulbourne [^1]

Parameters: Gc, Ge, λmax

Model: ``\\frac{G_c}{6}I_1-G_c\\lambda_{max}\\log\\bigg(3\\lambda_{max}^2-I_1\\bigg)+G_e\\sum\\limits_{i=1}^{3}\\big(\\lambda_i+\\frac{1}{\\lambda_i}\\big)``

[^1]: > Davidson JD, Goulbourne NC. A nonaffine network model for elastomers undergoing finite deformations. Journal of the Mechanics and Physics of Solids. 2013 Aug 1;61(8):1784-97.
"""
struct DavidsonGoulbourne <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::DavidsonGoulbourne, λ⃗::AbstractVector, (; Gc, Ge, λmax))
    1 / 6 * Gc * I₁(λ⃗) - Gc * λmax^2 * log(3*λmax^2 - I₁(λ⃗)) + Ge * (λ⃗[1] + 1 / λ⃗[1] + λ⃗[2] + 1 / λ⃗[2] + λ⃗[3] + 1 / λ⃗[3])
end

function parameters(ψ::DavidsonGoulbourne)
    return (:Gc, :Ge, :λmax)
end

function parameter_bounds(ψ::DavidsonGoulbourne, data::AbstractHyperelasticData)
    I₁_max = maximum(I₁.(collect.(data.λ⃗)))
    λmax_min = sqrt(I₁_max / 3)
    lb = (Gc=0, Ge=0, λmax=λmax_min)
    ub = nothing
    return (lb=lb, ub=ub)
end

"""
Khiêm-Itskov Model [^1]

Parameters: μcκ, n, q, μt

Model: ``\\mu_c \\kappa n \\log\\bigg(\\frac{\\sin(\\frac{\\pi}{\\sqrt{n}})(\\frac{I_1}{3})^{\\frac{q}{2}}}{\\sin(\\frac{\\pi}{\\sqrt{n}}(\\frac{I_1}{3})^{\\frac{q}{2}}}\\bigg)+\\mu_t\\big[\\frac{I_2}{3}^{1/2} - 1 \\big]``

[^1]: > Khiêm VN, Itskov M. Analytical network-averaging of the tube model:: Rubber elasticity. Journal of the Mechanics and Physics of Solids. 2016 Oct 1;95:254-69.
"""
struct KhiemItskov <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::KhiemItskov, λ⃗::AbstractVector, (; μcκ, n, q, μt))
    μcκ * n * log((sin(π / sqrt(n)) * (I₁(λ⃗) / 3)^(q / 2)) / (sin(π / sqrt(n) * (I₁(λ⃗) / 3)^(q / 2)))) + μt * ((I₂(λ⃗) / 3)^(1 / 2) - 1)
end

function ContinuumModels.StrainEnergyDensity(ψ::KhiemItskov, I⃗, (; μcκ, n, q, μt), I::InvariantForm)
    num = (sin(π / sqrt(n)) * (I⃗[1] / 3)^(q / 2))
    denom = (sin(π / sqrt(n) * (I⃗[1] / 3)^(q / 2)))
    @assert num ≥ denom "Parameters are not feasible"
    μcκ * n * log(num / denom) + μt * ((I⃗[2] / 3)^(1 / 2) - 1)
end

function parameters(ψ::KhiemItskov)
    return (:μcκ, :n, :q, :μt)
end


function constraints(ψ::KhiemItskov, data::AbstractHyperelasticData)
    I₁_max = maximum(I₁.(data.λ⃗))
    f(u, p) = [(sin(π / sqrt(u.n)) * (I₁_max / 3)^(u.q / 2)) / (sin(π / sqrt(u.n) * (I₁_max / 3)^(u.q / 2)))]
    return f
end


struct GeneralConstitutiveModel_Network <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::GeneralConstitutiveModel_Network, λ⃗::AbstractVector, (; Gc, N))
    I1 = I₁(λ⃗)
    Gc * N * log((3 * N + 0.5 * I1) / (3 * N - I1))
end

function parameters(ψ::GeneralConstitutiveModel_Network)
    return (:Gc, :N)
end

function parameter_bounds(ψ::GeneralConstitutiveModel_Network, data::AbstractHyperelasticData)
    I₁_max = maximum(I₁.(collect.(data.λ⃗)))
    N_min = I₁_max / 3
    lb = (Gc=-Inf, N=N_min)
    ub = nothing
    return (lb=lb, ub=ub)
end

struct GeneralConstitutiveModel_Tube <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::GeneralConstitutiveModel_Tube, λ⃗::AbstractVector, (; Ge))
    @tullio W := Ge / λ⃗[i]
end

function parameters(ψ::GeneralConstitutiveModel_Tube)
    return (:Ge,)
end

function parameter_bounds(ψ::GeneralConstitutiveModel_Tube, data::AbstractHyperelasticData)
    lb = nothing
    ub = nothing
    return (lb=lb, ub=ub)
end

"""
General Constitutive Model [^1]

Parameters: Gc, Ge, N

Model: ``G_c N \\log\\bigg(\\frac{3N+\\frac{1}{2}I_1}{3N-I_1}\\bigg)+G_e\\sum\\limits_{i=1}^{3}\\frac{1}{\\lambda_I}``

[^1]: > Xiang Y, Zhong D, Wang P, Mao G, Yu H, Qu S. A general constitutive model of soft elastomers. Journal of the Mechanics and Physics of Solids. 2018 Aug 1;117:110-22.
"""
struct GeneralConstitutiveModel <: AbstractHyperelasticModel end

function ContinuumModels.StrainEnergyDensity(ψ::GeneralConstitutiveModel, λ⃗::AbstractVector, ps)
    ContinuumModels.StrainEnergyDensity(GeneralConstitutiveModel_Network(), λ⃗, ps) + ContinuumModels.StrainEnergyDensity(GeneralConstitutiveModel_Tube(), λ⃗, ps)
end

function parameters(ψ::GeneralConstitutiveModel)
    return (:Gc, :Ge, :N)
end

function parameter_bounds(ψ::GeneralConstitutiveModel, data::AbstractHyperelasticData)
    I₁_max = maximum(I₁.(collect.(data.λ⃗)))
    N_min = I₁_max / 3
    lb = (Gc=-Inf, Ge=-Inf, N=N_min)
    ub = nothing
    return (lb=lb, ub=ub)
end


"""
Full Network - Wu Geisson [^1][^2][^3]

Parameters: μ, N, ρ

Model: ``(1-\\rho)W_{3Chain}+\\rho W_{8chain}``

[^1]: > Treloar LR, Riding G. A non-Gaussian theory for rubber in biaxial strain. I. Mechanical properties. Proceedings of the Royal Society of London. A. Mathematical and Physical Sciences. 1979 Dec 31;369(1737):261-80.
[^2]: > Wu PD, van der Giessen E. On improved 3-D non-Gaussian network models for rubber elasticity. Mechanics research communications. 1992 Sep 1;19(5):427-33.
[^3]: > Wu PD, Van Der Giessen E. On improved network models for rubber elasticity and their applications to orientation hardening in glassy polymers. Journal of the Mechanics and Physics of Solids. 1993 Mar 1;41(3):427-56.
"""
struct FullNetwork <: AbstractHyperelasticModel
    ℒinv::Function
    FullNetwork(; ℒinv::Function=TreloarApproximation) = new(ℒinv)
end

function ContinuumModels.StrainEnergyDensity(ψ::FullNetwork, λ⃗::AbstractVector, (; μ, N, ρ))
    W3 = ContinuumModels.StrainEnergyDensity(ThreeChainModel(ℒinv=ψ.ℒinv), λ⃗, (μ=μ, N=N))
    W8 = ContinuumModels.StrainEnergyDensity(ArrudaBoyce(ℒinv=ψ.ℒinv), λ⃗, (μ=μ, N=N))
    (1 - ρ) * W3 + ρ * W8
end

function parameters(ψ::FullNetwork)
    return (:μ, :N, :ρ)
end

function parameter_bounds(ψ::FullNetwork, data::AbstractHyperelasticData)
    I₁_max = maximum(I₁.(collect.(data.λ⃗)))
    λ_max = maximum(maximum.(collect.(data.λ⃗)))
    N₁ = λ_max^2
    N₂ = I₁_max / 3
    N_min = (N₁ > N₂) ? N₁ : N₂
    lb = (μ=-Inf, N=N_min, ρ=-Inf)
    ub = nothing
    return (lb=lb, ub=ub)
end

"""
Zuniga - Beatty [^1]

Parameters: μ, N₃, N₈

Model: ``\\sqrt{\\frac{N_3+N_8}{2N_3}}W_{3Chain}+\\sqrt{\\frac{I_1}{3N_8}}W_{8Chain}``

[^1]: > Elı́as-Zúñiga A, Beatty MF. Constitutive equations for amended non-Gaussian network models of rubber elasticity. International journal of engineering science. 2002 Dec 1;40(20):2265-94.
"""
struct ZunigaBeatty <: AbstractHyperelasticModel
    ℒinv::Function
    ZunigaBeatty(; ℒinv::Function=TreloarApproximation) = new(ℒinv)
end

function ContinuumModels.StrainEnergyDensity(ψ::ZunigaBeatty, λ⃗::AbstractVector, (; μ, N₃, N₈))
    ΛL = √((N₃ + N₈) / 2)
    ρ₃ = ΛL / √(N₃)
    W3 = ContinuumModels.StrainEnergyDensity(ThreeChainModel(ℒinv=ψ.ℒinv), λ⃗, (μ=μ, N=N₃))
    W8 = ContinuumModels.StrainEnergyDensity(ArrudaBoyce(ℒinv=ψ.ℒinv), λ⃗, (μ=μ, N=N₈))
    Λch = 1 / √(3) * √(I₁(λ⃗))
    ρ₈ = Λch / √(N₈)
    return ρ₃ * W3 + ρ₈ * W8
end

function parameters(ψ::ZunigaBeatty)
    return (:μ, :N₃, :N₈)
end

function parameter_bounds(ψ::ZunigaBeatty, data::AbstractHyperelasticData)
    λ_max = maximum(maximum.(collect.(data.λ⃗)))
    I₁_max = maximum(I₁.(collect.(data.λ⃗)))
    N₃_min = λ_max^2
    N₈_min = I₁_max / 3
    lb = (μ=-Inf, N₃=N₃_min, N₈=N₈_min)
    ub = nothing
    return (lb=lb, ub=ub)
end
"""
Lim [^1]

Parameters: μ₁, μ₂, N, Î₁

Model: ``(1-f(\\frac{I_1-3}{\\hat{I_1}-3}))W_{NeoHookean}(μ₁)+fW_{ArrudaBoyce}(μ₂, N)``

[^1]: > Lim GT. Scratch behavior of polymers. Texas A&M University; 2005.
"""
struct Lim <: AbstractHyperelasticModel
    ℒinv::Function
    Lim(; ℒinv::Function=TreloarApproximation) = new(ℒinv)
end

function ContinuumModels.StrainEnergyDensity(ψ::Lim, λ⃗::AbstractVector, (; μ₁, μ₂, N, Î₁))
    Wg = ContinuumModels.StrainEnergyDensity(NeoHookean(), λ⃗, (μ=μ₁,))
    W8 = ContinuumModels.StrainEnergyDensity(ArrudaBoyce(ℒinv=ψ.ℒinv), λ⃗, (μ=μ₂, N=N))
    f(x) = x^3 * (10 - 15x + 6x^2)
    ζ = (I₁(λ⃗) - 3) / (Î₁ - 3)
    (1 - f(ζ)) * Wg + f(ζ) * W8
end

function ContinuumModels.StrainEnergyDensity(ψ::Lim, I⃗, (; μ₁, μ₂, N, Î₁), I::InvariantForm)
    Wg = ContinuumModels.StrainEnergyDensity(NeoHookean(), I⃗, (μ = μ₁), I)
    W8 = ContinuumModels.StrainEnergyDensity(ArrudaBoyce(ℒinv=ψ.ℒinv), I⃗, (μ=μ₂, N=N), I)
    f(x) = x^3 * (10 - 15x + 6x^2)
    ζ = (I⃗[1] - 3) / (Î₁ - 3)
    (1 - f(ζ)) * Wg + f(ζ) * W8
end

function parameters(ψ::Lim)
    return (:μ₁, :μ₂, :N, :Î₁)
end

function parameter_bounds(ψ::Lim, data::AbstractHyperelasticData)
    I₁_max = maximum(I₁.(collect.(data.λ⃗)))
    N_min = I₁_max / 3
    lb = (μ₁=-Inf, μ₂=-Inf, N=N_min, Î₁=3)
    ub = nothing
    return (lb=lb, ub=ub)
end

"""
Bechir Chevalier [^1]

Parameters: μ₀, η, ρ, N₃, N₈

Model:

``W_{3Chain}(\\mu_f, N_3)+W_{8Chain}(\\frac{\\mu_c}{3}, N_8)``

``\\mu_f = \\rho\\sqrt{\\frac{I_1}{3N_8}}``

``\\mu_c = \\bigg(1-\\frac{\\eta\\alpha}{\\sqrt{N_3}}\\bigg)\\mu_0``

``\\alpha = \\max{\\lambda_1, \\lambda_2, \\lambda_3}``

[^1]: > Bechir H, Chevalier L, Idjeri M. A three-dimensional network model for rubber elasticity: The effect of local entanglements constraints. International journal of engineering science. 2010 Mar 1;48(3):265-74.
"""
struct BechirChevalier <: AbstractHyperelasticModel
    ℒinv::Function
    BechirChevalier(; ℒinv::Function=TreloarApproximation) = new(ℒinv)
end

function ContinuumModels.StrainEnergyDensity(ψ::BechirChevalier, λ⃗::AbstractVector, (; μ₀, η, ρ, N₃, N₈))
    μf = ρ * √(I₁(λ⃗) / 3 / N₈)
    W3 = ContinuumModels.StrainEnergyDensity(ThreeChainModel(ℒinv=ψ.ℒinv), λ⃗, (μ=μf, N=N₃))
    α = maximum(λ⃗)
    μc = (1 - η * α / √(N₃)) * μ₀
    W8 = ContinuumModels.StrainEnergyDensity(ArrudaBoyce(ℒinv=ψ.ℒinv), λ⃗, (μ=μc / 3, N=N₈))
    W3 + W8
end

function parameters(ψ::BechirChevalier)
    return (:μ₀, :η, :ρ, :N₃, :N₈)
end

function parameter_bounds(ψ::BechirChevalier, data::AbstractHyperelasticData)
    lb = (μ₀=-Inf, η=-Inf, ρ=-Inf, N₃=0, N₈=0)
    ub = nothing
    return (lb=lb, ub=ub)
end
