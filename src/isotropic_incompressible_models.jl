# # Available Models
export MooneyRivlin,
    NeoHookean,
    Gent,
    Biderman,
    Isihara,
    JamesGreenSimpson,
    Lion,
    Yeoh,
    HauptSedlan,
    HartmannNeff,
    HainesWilson,
    Carroll,
    BahremanDarijani,
    Zhao,
    Knowles,
    Swanson,
    YamashitaKawabata,
    DavisDeThomas,
    Gregory,
    ModifiedGregory,
    Beda,
    Amin,
    LopezPamies,
    GenYeoh,
    VerondaWestmann,
    FungDemiray,
    Vito,
    ModifiedYeoh,
    MansouriDarijani,
    GentThomas,
    HossMarczakI,
    HossMarczakII,
    ExpLn,
    VanDerWaals,
    TakamizawaHayashi,
    YeohFleming,
    PucciSaccomandi,
    HorganSaccomandi,
    Beatty,
    ArrudaBoyce,
    Ogden,
    EdwardVilgis,
    NonaffineTube,
    Tube,
    MCC,
    Bechir4Term,
    ConstrainedJunction,
    ContinuumHybrid,
    ArmanNarooei,
    PengLandel,
    ValanisLandel,
    Attard,
    Shariff,
    ThreeChainModel,
    ModifiedFloryErman,
    ABGI,
    BechirChevalier,
    Bootstrapped8Chain,
    DavidsonGoulbourne,
    ExtendedTubeModel,
    FullNetwork,
    HartSmith,
    GeneralConstitutiveModel,
    Lim,
    NonaffineMicroSphere,
    AffineMicroSphere,
    ZunigaBeatty,
    ChevalierMarco,
    Alexander,
    GornetDesmorat,
    LambertDianiRey,
    AnsarriBenam

export HorganMurphy, KhiemItskov

export GeneralConstitutiveModel_Network, GeneralConstitutiveModel_Tube

struct ArrudaBoyce{T} <: AbstractIncompressibleModel{T}
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation
end

"""
$(SIGNATURES)

# Model:

```math
W = \\mu N \\left( \\frac{\\lambda_{chain}}{\\sqrt{N}} \\beta + \\log\\left(\\frac{\\beta}{\\sinh\\beta}\\right) \\right)
```

where

```math
\\beta = \\mathcal{L}^{-1}\\left(\\frac{\\lambda_{chain}}{\\sqrt{N}}\\right)
```

and

```math
\\lambda_{chain} = \\sqrt{\\frac{I_1}{3}}
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()`
- `ℒinv=TreloarApproximation()`: Sets the inverse Langevin approxamation used

# Parameters:
- `μ`: Small strain shear modulus
- `N`: Square of the locking stretch of the network.

> Arruda EM, Boyce MC. A three-dimensional constitutive model for the large stretch behavior of rubber elastic materials. Journal of the Mechanics and Physics of Solids. 1993 Feb 1;41(2):389-412.

"""
ArrudaBoyce(
    type::T = PrincipalValueForm();
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation = TreloarApproximation(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = ArrudaBoyce{T}(ℒinv)

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::ArrudaBoyce,
    λ⃗::Vector{T},
    (; μ, N),
) where {T}
    rchain_Nl = √(I₁(λ⃗) / 3 / N)
    β = inverse_langevin_approximation(rchain_Nl, ψ.ℒinv)
    return μ * N * (rchain_Nl * β + log(β / sinh(β)))
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::ArrudaBoyce{I},
    I⃗::Vector{T},
    (; μ, N),
) where {T,I<:InvariantForm}
    rchain_Nl = √(I⃗[1] / 3 / N)
    β = inverse_langevin_approximation(rchain_Nl, ψ.ℒinv,)
    return μ * N * (rchain_Nl * β + log(β / sinh(β)))
end

function parameters(::ArrudaBoyce)
    return (:μ, :N)
end

function parameter_bounds(::ArrudaBoyce, data::AbstractHyperelasticTest)
    I₁_max = maximum(I₁.(data.data.λ))
    N_max = I₁_max / 3
    lb = (μ = -Inf, N = N_max)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct ABGI{T} <: AbstractIncompressibleModel{T}
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation
    AB::ArrudaBoyce
end

"""
$(SIGNATURES)

# Model:

```math
W = W_{Arruda-Boyce} + \\frac{G_e}{n} \\left(\\sum_{i=1}^{3}\\lambda_i^n-3\\right)
```

# Arguments:

- `ℒinv = TreloarApproximation()`: Sets the inverse Langevin approxamationused (default = `TreloarApproximation()`)

# Parameters:

- `μ`
- `N`
- `Ge`
- `n`

> Meissner B, Matějka L. A Langevin-elasticity-theory-based constitutiveequation for rubberlike networks and its comparison with biaxialstress–strain data. Part I. Polymer. 2003 Jul 1;44(16):4599-610.
"""
ABGI(;
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation = TreloarApproximation(),
) = ABGI{PrincipalValueForm}(ℒinv, ArrudaBoyce(PrincipalValueForm(), ℒinv = ℒinv))

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::ABGI,
    λ⃗::Vector{T},
    (; μ, N, Ge, n),
) where {T}
    WAB = StrainEnergyDensity(ψ.AB, λ⃗, (μ = μ, N = N))
    WGI = Ge * (λ⃗[1]^n + λ⃗[2]^n + λ⃗[3]^n - 3) / n
    return WAB + WGI
end

function parameters(::ABGI)
    return (:μ, :N, :Ge, :n)
end

function parameter_bounds(::ABGI, data::AbstractHyperelasticTest)
    I₁_max = maximum(I₁.(data.data.λ))
    lb = (μ = -Inf, N = 11 / 35 * I₁_max, Ge = -Inf, n = 0.0)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct AffineMicroSphere{T,R,S} <: AbstractIncompressibleModel{T}
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation
    r⃗::Vector{R}
    w::Vector{S}
    λr::Function
end

"""
$(SIGNATURES)

# Model:
- See Paper

# Arguments:
- `ℒinv=TreloarApproximation()`: Sets the inverse Langevin approxamation used
- `n::Int = 21`: Number of quadrature points for the spherical integration

# Parameters:
- `μ`
- `N`

> Miehe C, Göktepe S, Lulei F. A micro-macro approach to rubber-like materials—part I: the non-affine micro-sphere model of rubber elasticity. Journal of the Mechanics and Physics of Solids. 2004 Nov 1;52(11):2617-60.
"""
function AffineMicroSphere(;
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation = TreloarApproximation(),
    n = 21,
)
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
    λr((; λ, N), r) = sqrt(sum(λ .^ 2 .* r .^ 2)) / √N
    AffineMicroSphere{PrincipalValueForm,eltype(r⃗),eltype(w)}(ℒinv, r⃗, w, λr)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::AffineMicroSphere,
    λ⃗::Vector{T},
    (; μ, N),
) where {T}
    λr = map(Base.Fix1(ψ.λr, (λ = λ⃗, N = N)), ψ.r⃗)
    β = @. inverse_langevin_approximation(λr, ψ.ℒinv)
    ψf = @. μ * N * (λr * β + log(β / sinh(β))) * ψ.w
    return sum(ψf)
end

function parameters(::AffineMicroSphere)
    return (:μ, :N)
end

function parameter_bounds(ψ::AffineMicroSphere, test::AbstractHyperelasticTest)
    λ = test.data.λ
    λr = maximum(x -> map(Base.Fix1(ψ.λr, (λ = x, N = 1)), ψ.r⃗), λ)
    N_min = maximum(λr)
    lb = (μ = -Inf, N = N_min)
    ub = nothing
    (lb = lb, ub = ub)
end

struct Alexander{T} <: AbstractIncompressibleModel{T} end

"""
$(SIGNATURES)

# Model:

```math
W = \\frac{C_1 \\sqrt{\\pi}\\text{erfi}\\big(\\sqrt{k}(I_1-3)\\big)}{2\\sqrt{k}}+C_2\\log{\\frac{I_2-3+\\gamma}{\\gamma}}+C_3(I_2-3)
```

# Parameters:
- `μ`
- `C₁`
- `C₂`
- `C₃`
- `k`
- `γ`

> Alexander H. A constitutive relation for rubber-like materials. International Journal of Engineering Science. 1968 Sep 1;6(9):549-63.
"""
Alexander() = Alexander{PrincipalValueForm}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Alexander{T},
    λ⃗::Vector{S},
    (; μ, C₁, C₂, C₃, k, γ),
) where {T<:PrincipalValueForm,S}
    return μ / 3 * (
        C₁ * √π * erfi(√k * (I₁(λ⃗) - 3)) / 2 / √k +
        C₂ * log((I₂(λ⃗) - 3 + γ) / γ) +
        C₃ * (I₂(λ⃗) - 3)
    )
end

function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    ::Alexander{T},
    λ⃗::Vector{S},
    (; μ, C₁, C₂, C₃, k, γ);
    kwargs...,
) where {T<:PrincipalValueForm,S}
    I1 = I₁(λ⃗)
    I2 = I₂(λ⃗)
    s = @. μ / 3 * (
        (3 * λ⃗^2 - I1) * C₁ * exp(k * (I1 - 3)^2) +
        (I2 - 3 * λ⃗^2) * (C₂ / (I2 - 3 + γ) + C₃)
    )
    return s
end

parameters(::Alexander) = (:μ, :C₁, :C₂, :C₃, :k, :γ)

struct MooneyRivlin{T} <: AbstractIncompressibleModel{T}
    GMR::GeneralMooneyRivlin{T}

end

"""
$(SIGNATURES)

# Model:

```math
W = C_{10}(I_1-3)+C_{01}(I_2-3)
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()`

# Parameters:
- `C01`
- `C10`

> Mooney M. A theory of large elastic deformation. Journal of applied physics. 1940 Sep;11(9):582-92.
"""
MooneyRivlin(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    MooneyRivlin{T}(GeneralMooneyRivlin(T()))

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::MooneyRivlin{T},
    I⃗::Vector{S},
    (; C10, C01),
) where {T,S}
    ContinuumMechanicsBase.StrainEnergyDensity(ψ.GMR, I⃗, (C⃗ = [
        0.0 C10
        C01 0.0
    ],))
end

parameters(::MooneyRivlin) = (:C10, :C01)

struct NeoHookean{T} <: AbstractIncompressibleModel{T} end

"""
$(SIGNATURES)

# Model:

```math
W = \\frac{\\mu}{2}(I_1-3)
```

# Arguments
- `type = PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `μ`: Small strain shear modulus

> Treloar LR. The elasticity of a network of long-chain molecules—II. Transactions of the Faraday Society. 1943;39:241-6.
"""
NeoHookean(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = NeoHookean{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::NeoHookean{T},
    λ⃗::Vector{S},
    (; μ),
) where {T<:PrincipalValueForm,S}
    μ / 2 * (I₁(λ⃗) - 3)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::NeoHookean{T},
    I⃗::Vector{S},
    (; μ),
) where {T<:InvariantForm,S}
    μ / 2 * (I⃗[1] - 3)
end

parameters(::NeoHookean) = (:μ,)

struct Isihara{T} <: AbstractIncompressibleModel{T}
    GMR::GeneralMooneyRivlin{T}
end

"""
$(SIGNATURES)

# Model:

```math
W = \\sum\\limits_{i,j=0}^{2, 1}C_{i,j}(I_1-3)^i(I_2-3)^j
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()`

# Parameters:
- `C10`
- `C20`
- `C01`

> Isihara A, Hashitsume N, Tatibana M. Statistical theory of rubber‐like elasticity. IV.(two‐dimensional stretching). The Journal of Chemical Physics. 1951 Dec;19(12):1508-12.
"""
Isihara(type::T = PrincipalValueForm()) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    Isihara{T}(GeneralMooneyRivlin(type))

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::Isihara{T},
    λ⃗::Vector{S},
    (; C10, C20, C01),
) where {T,S}
    StrainEnergyDensity(ψ.GMR, λ⃗, (C⃗ = [
        0.0 C10 C20
        C01 0.0 0.0
    ],))
end

parameters(ψ::Isihara) = (:C10, :C20, :C01)

struct Biderman{T} <: AbstractIncompressibleModel{T} end


"""
$(SIGNATURES)

# Model:

```math
W = \\sum\\limits_{i,j=0}^{3, 1}C_{i,j}(I_1-3)^i(I_2-3)^j
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()`

# Parameters:
- `C10`
- `C01`
- `C20`
- `C30`

> Biderman VL. Calculation of rubber parts. Rascheti na prochnost. 1958;40.
"""
Biderman(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = Biderman{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::Biderman{T},
    λ⃗::Vector{S},
    (; C10, C01, C20, C30),
) where {T<:PrincipalValueForm,S}
    I1 = I₁(λ⃗)
    I2 = I₂(λ⃗)
    W = C10 * (I1 - 3) + C20 * (I1 - 3)^2 + C30 * (I1 - 3)^3 + C01 * (I2 - 3)
    return W
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::Biderman{T},
    I⃗::Vector{S},
    (; C10, C01, C20, C30),
) where {T<:InvariantForm,S}
    I1 = I⃗[1]
    I2 = I⃗[2]
    W = C10 * (I1 - 3) + C20 * (I1 - 3)^2 + C30 * (I1 - 3)^3 + C01 * (I2 - 3)
    return W
end

parameters(::Biderman) = (:C10, :C01, :C20, :C30)

struct JamesGreenSimpson{T} <: AbstractIncompressibleModel{T}
    GMR::GeneralMooneyRivlin{T}
end


"""
$(SIGNATURES)

# Model:

```math
W = \\sum\\limits_{i,j=0}^{3, 1}C_{i,j}(I_1-3)^i(I_2-3)^j
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()`

# Parameters:
- `C10`
- `C01`
- `C11`
- `C20`
- `C30`

> James AG, Green A, Simpson GM. Strain energy functions of rubber. I. Characterization of gum vulcanizates. Journal of Applied Polymer Science. 1975 Jul;19(7):2033-58.
"""
JamesGreenSimpson(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    JamesGreenSimpson{T}(GeneralMooneyRivlin(type))

function ContinuumMechanicsBase.StrainEnergyDensity(
    W::JamesGreenSimpson{T},
    λ⃗::Vector{S},
    (; C10, C01, C11, C20, C30),
) where {T,S}
    StrainEnergyDensity(W.GMR, λ⃗, (C⃗ = [
        0.0 C10 C20 C30
        C01 0.0 0.0 0.0
    ],))
end

parameters(ψ::JamesGreenSimpson) = (:C10, :C01, :C11, :C20, :C30)

struct HainesWilson{T} <: AbstractIncompressibleModel{T}
    GMR::GeneralMooneyRivlin{T}
end

"""
$(SIGNATURES)

# Model:

```math
W = \\sum\\limits_{i,j=0}^{3, 2}C_{i,j}(I_1-3)^i(I_2-3)^j
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()`

# Parameters:
- `C10`
- `C01`
- `C11`
- `C02`
- `C20`
- `C30`

> Haines DW, Wilson WD. Strain-energy density function for rubberlike materials. Journal of the Mechanics and Physics of Solids. 1979 Aug 1;27(4):345-60.
"""
HainesWilson(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    HainesWilson{T}(GeneralMooneyRivlin(type))

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::HainesWilson{T},
    λ⃗::Vector{S},
    (; C10, C01, C11, C02, C20, C30),
) where {T,S}
    StrainEnergyDensity(ψ.GMR, λ⃗, (C⃗ = [
        0.0 C10 C20 C30
        C01 C11 0.0 0.0
        C02 0.0 0.0 0.0
    ],))
end

parameters(::HainesWilson) = (:C10, :C01, :C11, :C02, :C20, :C30)

struct Yeoh{T} <: AbstractIncompressibleModel{T}
    GMR::GeneralMooneyRivlin{T}
end

"""
$(SIGNATURES)

# Model:

```math
W = \\sum\\limits_{i,j=0}^{3, 0}C_{i,j}(I_1-3)^i(I_2-3)^j
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()`

# Parameters:
- C10
- C20
- C30

> Yeoh OH. Characterization of elastic properties of carbon-black-filled rubber vulcanizates. Rubber chemistry and technology. 1990 Nov;63(5):792-805.
"""
Yeoh(type::T = PrincipalValueForm()) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    Yeoh{T}(GeneralMooneyRivlin(type))

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::Yeoh{T},
    λ⃗::Vector{S},
    (; C10, C20, C30),
) where {T,S}
    StrainEnergyDensity(ψ.GMR, λ⃗, (C⃗ = [0.0 C10 C20 C30],))
end

parameters(ψ::Yeoh) = (:C10, :C20, :C30)

struct Lion{T} <: AbstractIncompressibleModel{T}
    GMR::GeneralMooneyRivlin{T}
end

"""
$(SIGNATURES)

# Model:

```math
W = \\sum\\limits_{i,j=0}^{5,1}C_{i,j}(I_1-3)^i(I_2-3)^j
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()`

# Parameters:
- `C10`
- `C01`
- `C50`

> Lion A. On the large deformation behaviour of reinforced rubber at different temperatures. Journal of the Mechanics and Physics of Solids. 1997 Nov 1;45(11-12):1805-34.
"""
Lion(type::T = PrincipalValueForm()) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    Lion{T}(GeneralMooneyRivlin(type))

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::Lion{T},
    λ⃗::Vector{S},
    (; C10, C01, C50),
) where {T,S}
    StrainEnergyDensity(ψ.GMR, λ⃗, (C⃗ = [
        0.0 C10 0.0 0.0 0.0 C50
        C01 0.0 0.0 0.0 0.0 0.0
    ],))
end

parameters(::Lion) = (:C10, :C01, :C50)

struct HauptSedlan{T} <: AbstractIncompressibleModel{T}
    GMR::GeneralMooneyRivlin{T}
end

"""
$(SIGNATURES)

# Model:

```math
W = \\sum\\limits_{i,j=0}^{3, 2}C_{i,j}(I_1-3)^i(I_2-3)^j
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()`

# Parameters:
- `C10`
- `C01`
- `C11`
- `C02`
- `C30`

> Haupt P, Sedlan K. Viscoplasticity of elastomeric materials: experimental facts and constitutive modelling. Archive of Applied Mechanics. 2001 Mar;71(2):89-109.
"""
HauptSedlan(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    HauptSedlan{T}(GeneralMooneyRivlin(type))

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::HauptSedlan{T},
    λ⃗::Vector{S},
    (; C10, C01, C11, C02, C30),
) where {T,S}
    StrainEnergyDensity(ψ.GMR, λ⃗, (C⃗ = [
        0.0 C10 0.0 C30
        C01 C11 0.0 0.0
        C02 0.0 0.0 0.0
    ],))
end

parameters(::HauptSedlan) = (:C10, :C01, :C11, :C02, :C30)

struct HartmannNeff{T} <: AbstractIncompressibleModel{T} end

"""
$(SIGNATURES)

# Model:

```math
W = \\sum\\limits_{i,j=0}^{M,N}C_{i,0}(I_1-3)^i -3\\sqrt{3}^j+\\alpha(I_1-3)
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()`

# Parameters:
- `α`
- `Ci⃗0`
- `C0j⃗`

> Hartmann S, Neff P. Polyconvexity of generalized polynomial-type hyperelastic strain energy functions for near-incompressibility. International journal of solids and structures. 2003 Jun 1;40(11):2767-91.
"""
HartmannNeff(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = HartmannNeff{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::HartmannNeff{T},
    λ⃗::Vector{S},
    (; α, Ci⃗0, C0j⃗),
) where {T<:PrincipalValueForm,S}
    I1 = I₁(λ⃗)
    I2 = I₂(λ⃗)
    i_max = length(Ci⃗0)
    j_max = length(C0j⃗)
    W1 = @. Ci⃗0 * (I1 - 3)^(1:i_max)
    W2 = @. C0j⃗ * (I2^(3 / 2) - 3sqrt(3))^(1:j_max)
    return sum(W1) + sum(W2) + α * (I1^3 - 3^3)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::HartmannNeff{T},
    I⃗::Vector{S},
    (; α, Ci⃗0, C0j⃗),
) where {T<:InvariantForm,S}
    i_max = length(Ci⃗0)
    j_max = length(C0j⃗)
    W1 = @. Ci⃗0 * (I⃗[1] - 3)^(1:i_max)
    W2 = @. C0j⃗ * (I⃗[2]^(3 / 2) - 3sqrt(3))^(1:j_max)
    return sum(W1) + sum(W2) + α * (I⃗[1]^3 - 3^3)
end

parameters(::HartmannNeff) = (:α, :Ci⃗0, :C0j⃗)

struct Carroll{T} <: AbstractIncompressibleModel{T} end

"""
$(SIGNATURES)

# Model:

```math
W = AI_1+BI_1^4+C\\sqrt{I_2}
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()`

# Parameters:
- `A`
- `B`
- `C`

> Carroll M. A strain energy function for vulcanized rubbers. Journal of Elasticity. 2011 Apr;103(2):173-87.
"""
Carroll(type::T = PrincipalValueForm()) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    Carroll{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Carroll{T},
    λ⃗::Vector{S},
    (; A, B, C),
) where {T<:PrincipalValueForm,S}
    return A * I₁(λ⃗) + B * I₁(λ⃗)^4 + C * I₂(λ⃗)^(1 / 2)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Carroll{T},
    I⃗::Vector{S},
    (; A, B, C),
) where {T<:InvariantForm,S}
    return A * I⃗[1] + B * I⃗[1]^4 + C * I⃗[2]^(1 / 2)
end

parameters(::Carroll) = (:A, :B, :C)

struct BahremanDarijani{PrincipalValueForm} <:
       AbstractIncompressibleModel{PrincipalValueForm}
    GDN::GeneralDarijaniNaghdabadi
end

"""
$(SIGNATURES)

# Model:

```math
W = \\sum\\limits_{i = 1}{3}\\sum\\limits_{j=0}^{N} A_j (\\lambda_i^{m_j}-1) + B_j(\\lambda_i^{-n_j}-1)
```
# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()`

# Parameters:
- `A2`
- `B2`
- `A4`
- `A6`

> Bahreman M, Darijani H. New polynomial strain energy function; application to rubbery circular cylinders under finite extension and torsion. Journal of Applied Polymer Science. 2015 Apr 5;132(13).
"""
BahremanDarijani(type::T = PrincipalValueForm()) where {T<:PrincipalValueForm} =
    BahremanDarijani{PrincipalValueForm}(GeneralDarijaniNaghdabadi(type))

function ContinuumMechanicsBase.StrainEnergyDensity(
    W::BahremanDarijani{T},
    λ⃗::Vector{S},
    (; A2, B2, A4, A6),
) where {T,S}
    StrainEnergyDensity(
        W.GDN,
        λ⃗,
        (A⃗ = [0, A2, 0, A4, 0, A6], B⃗ = [0, B2], m⃗ = [0, 2, 0, 4, 0, 6], n⃗ = [0, 2]),
    )
end

parameters(::BahremanDarijani) = (:A2, :B2, :A4, :A6)

struct Zhao{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = C_{-1}^{1}*(I_2-3)+C_{1}^{1}(I_1-3)+C_{2}^{1}(I_1^2-2I_2-3)+C_{2}^{2}(I_1^2-2I_2-3)^2
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()`

# Parameters:
- `C₋₁¹`
- `C₁¹`
- `C₂¹`
- `C₂²`
> Zhao Z, Mu X, Du F. Modeling and verification of a new hyperelastic modelfor rubber-like materials. Mathematical Problems in Engineering. 2019 May 22019.
"""
Zhao(type::T = PrincipalValueForm()) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    Zhao{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Zhao{T},
    λ⃗::Vector{S},
    (; C₋₁¹, C₁¹, C₂¹, C₂²),
) where {T<:PrincipalValueForm,S}
    return C₋₁¹ * (I₂(λ⃗) - 3) +
           C₁¹ * (I₁(λ⃗) - 3) +
           C₂¹ * (I₁(λ⃗)^2 - 2I₂(λ⃗) - 3) +
           C₂² * (I₁(λ⃗)^2 - 2I₂(λ⃗) - 3)^2
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Zhao{T},
    I⃗::Vector{S},
    (; C₋₁¹, C₁¹, C₂¹, C₂²),
) where {T<:InvariantForm,S}
    return C₋₁¹ * (I⃗[2] - 3) +
           C₁¹ * (I⃗[1] - 3) +
           C₂¹ * (I⃗[1]^2 - 2I⃗[2] - 3) +
           C₂² * (I⃗[1]^2 - 2I⃗[2] - 3)^2
end

parameters(::Zhao) = (:C₋₁¹, :C₁¹, :C₂¹, :C₂²)

struct Knowles{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = \\frac{\\mu}{2b}\\left(\\left(1+\\frac{b}{n}(I_1-3)\\right)^n-1\\right)
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- μ
- b
- n

> Knowles JK. The finite anti-plane shear field near the tip of a crack for a class of incompressible elastic solids. International Journal of Fracture. 1977 Oct;13(5):611-39.
"""
Knowles(type::T = PrincipalValueForm()) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    Knowles{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Knowles{T},
    λ⃗::Vector{S},
    (; μ, b, n),
) where {T<:PrincipalValueForm,S}
    return μ / (2b) * ((1 + (b / n) * (I₁(λ⃗) - 3))^n - 1)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Knowles{T},
    I⃗::Vector{S},
    (; μ, b, n),
) where {T<:InvariantForm,S}
    return μ / (2b) * ((1 + (b / n) * (I⃗[1] - 3))^n - 1)
end

parameters(::Knowles) = (:μ, :b, :n)

function parameter_bounds(::Knowles, data::AbstractHyperelasticTest)
    lb = (μ = -Inf, b = 0, n = 0)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct Swanson{T} <: AbstractIncompressibleModel{T} end

"""
$(SIGNATURES)

# Model:

```math
W = \\sum\\limits_{i=1}^{N} \\frac{3}{2}(\\frac{A_i}{1+\\alpha_i}(\\frac{I_1}{3})^{1+\\alpha_i}+\\frac{B_i}{1+\\beta_i}(\\frac{I_2}{3})^{1+\\beta_i}
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `A⃗`
- `α⃗`
- `B⃗`
- `β⃗`

> Swanson SR. A constitutive model for high elongation elastic materials.
"""
Swanson(type::T = PrincipalValueForm()) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    Swanson{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Swanson{T},
    λ⃗::Vector{S},
    (; A⃗, α⃗, B⃗, β⃗),
) where {T<:PrincipalValueForm,S}
    @assert length(A⃗) == length(α⃗) == length(B⃗) == length(β⃗) "The vectors are not the same length"
    I1 = I₁(λ⃗)
    I2 = I₂(λ⃗)
    return sum(
        @. 3 / 2 * (A⃗ / (1 + α⃗) * (I1 / 3)^(1 + α⃗) + B⃗ / (1 + β⃗) * (I2 / 3)^(1 + β⃗))
    )
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Swanson{T},
    I⃗::Vector{S},
    (; A⃗, α⃗, B⃗, β⃗),
) where {T<:InvariantForm,S}
    @assert length(A⃗) == length(α⃗) == length(B⃗) == length(β⃗) "The vectors are not the same length"
    return sum(
        @. 3 / 2 *
           (A⃗ / (1 + α⃗) * (I⃗[1] / 3)^(1 + α⃗) + B⃗ / (1 + β⃗) * (I⃗[2] / 3)^(1 + β⃗))
    )
end

parameters(::Swanson) = (:A⃗, :α⃗, :B⃗, :β⃗)

struct YamashitaKawabata{T} <: AbstractIncompressibleModel{T} end

"""
$(SIGNATURES)

# Model:

```math
W = C_1(I_1-3)+C_2(I_2-3)+\\frac{C_3}{N+1}(I_1-3)^{N+1}
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- C1
- C2
- C3
- N

> Yamashita Y, Kawabata S. Approximated form of the strain energy-density function of carbon-black filled rubbers for industrial applications. Nippon Gomu Kyokaishi(Journal of the Society of Rubber Industry, Japan)(Japan). 1992;65(9):517-28.
"""
YamashitaKawabata(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = YamashitaKawabata{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::YamashitaKawabata{T},
    λ⃗::Vector{S},
    (; C1, C2, C3, N),
) where {T<:PrincipalValueForm,S}
    return C1 * (I₁(λ⃗) - 3) + C2 * (I₂(λ⃗) - 3) + C3 / (N + 1) * (I₁(λ⃗) - 3)^(N + 1)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::YamashitaKawabata{T},
    I⃗::Vector{S},
    (; C1, C2, C3, N),
) where {T<:InvariantForm,S}
    return C1 * (I⃗[1] - 3) + C2 * (I⃗[2] - 3) + C3 / (N + 1) * (I⃗[1] - 3)^(N + 1)
end

parameters(::YamashitaKawabata) = (:C1, :C2, :C3, :N)


struct DavisDeThomas{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = \\frac{A}{2(1-\\frac{n}{2})}(I_1-3+C^2)^{1-\\frac{n}{2}}+k(I_1-3)^2
```

# Arguments
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `A`
- `n`
- `C`
- `k`

> Davies CK, De DK, Thomas AG. Characterization of the behavior of rubber for engineering design purposes. 1. Stress-strain relations. Rubber chemistry and technology. 1994 Sep;67(4):716-28.
"""
DavisDeThomas(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = DavisDeThomas{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::DavisDeThomas{T},
    λ⃗::Vector{S},
    (; A, n, C, k),
) where {T<:PrincipalValueForm,S}
    return A / (2 * (1 - n / 2)) * (I₁(λ⃗) - 3 + C^2)^(1 - n / 2) + k * (I₁(λ⃗) - 3)^2
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::DavisDeThomas{T},
    I⃗::Vector{S},
    (; A, n, C, k),
) where {T<:InvariantForm,S}
    return A / (2 * (1 - n / 2)) * (I⃗[1] - 3 + C^2)^(1 - n / 2) + k * (I⃗[1] - 3)^2
end

function parameters(::DavisDeThomas)
    return (:A, :n, :C, :k)
end

struct Gregory{T} <: AbstractIncompressibleModel{T} end

"""
$(SIGNATURES)

# Model:

```math
W = \\frac{A}{2-n}(I_1-3+C^2)^{1-\\frac{n}{2}}+\\frac{B}{2+m}(I_1-3+C^2)^{1+\\frac{m}{2}}
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `A`
- `B`
- `C`
- `m`
- `n`

> Gregory IH, Muhr AH, Stephens IJ. Engineering applications of rubber in simple extension. Plastics rubber and composites processing and applications. 1997;26(3):118-22.
"""
Gregory(type::T = PrincipalValueForm()) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    Gregory{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Gregory{T},
    λ⃗::Vector{S},
    (; A, B, C, m, n),
) where {T<:PrincipalValueForm,S}
    return A / (2 - n) * (I₁(λ⃗) - 3 + C^2)^(1 - n / 2) +
           B / (2 + m) * (I₁(λ⃗) - 3 + C^2)^(1 + m / 2)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Gregory{T},
    I⃗::Vector{S},
    (; A, B, C, m, n),
) where {T<:InvariantForm,S}
    return A / (2 - n) * (I⃗[1] - 3 + C^2)^(1 - n / 2) +
           B / (2 + m) * (I⃗[1] - 3 + C^2)^(1 + m / 2)
end

function parameters(::Gregory)
    return (:A, :B, :C, :m, :n)
end

struct ModifiedGregory{T} <: AbstractIncompressibleModel{T} end

"""
$(SIGNATURES)

# Model:

```math
W = \\frac{A}{1+\\alpha}(I_1-3+M^2)^{1+\\alpha}+\\frac{B}{1+\\beta}(I_1-3+N^2)^{1+\\beta}
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `A`
- `α`
- `M`
- `B`
- `β`
- `N`

> He H, Zhang Q, Zhang Y, Chen J, Zhang L, Li F. A comparative study of 85 hyperelastic constitutive models for both unfilled rubber and highly filled rubber nanocomposite material. Nano Materials Science. 2021 Jul 16.
"""
ModifiedGregory(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = ModifiedGregory{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::ModifiedGregory{T},
    λ⃗::Vector{S},
    (; A, α, M, B, β, N),
) where {T<:PrincipalValueForm,S}
    return A / (1 + α) * (I₁(λ⃗) - 3 + M^2)^(1 + α) +
           B / (1 + β) * (I₁(λ⃗) - 3 + N^2)^(1 + β)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::ModifiedGregory{T},
    I⃗::Vector{S},
    (; A, α, M, B, β, N),
) where {T<:InvariantForm,S}
    return A / (1 + α) * (I⃗[1] - 3 + M^2)^(1 + α) + B / (1 + β) * (I⃗[1] - 3 + N^2)^(1 + β)
end

function parameters(::ModifiedGregory)
    return (:A, :α, :M, :B, :β, :N)
end

struct Beda{T} <: AbstractIncompressibleModel{T}
    GB::GeneralBeda
end

"""
$(SIGNATURES)

# Model:

```math
W = \\frac{C_1}{\\alpha}(I_1-3)^{\\alpha}+C_2(I_1-3)+\\frac{C_3}{\\zeta}(I_1-3)^{\\zeta}+\\frac{K_1}{\\beta}(I_2-3)^\\beta
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `C1`
- `C2`
- `C3`
- `K1`
- `α`
- `β`
- `ζ`

> Beda T. Reconciling the fundamental phenomenological expression of the strain energy of rubber with established experimental facts. Journal of Polymer Science Part B: Polymer Physics. 2005 Jan 15;43(2):125-34.
"""
Beda(type::T = PrincipalValueForm()) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    Beda{T}(GeneralBeda(type))

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::Beda{T},
    λ⃗::Vector{S},
    (; C1, C2, C3, K1, α, β, ζ),
) where {T<:PrincipalValueForm,S}
    I1 = I₁(λ⃗)
    I2 = I₂(λ⃗)
    W1 = C1 / α * (I1 - 3)^α + C2 * (I1 - 3) + C3 / ζ * (I1 - 3)^ζ
    W2 = K1 / β * (I2 - 3)^β
    W = W1 + W2
    return W
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::Beda{T},
    I⃗::Vector{S},
    (; C1, C2, C3, K1, α, β, ζ),
) where {T<:InvariantForm,S}
    I1 = I⃗[1]
    I2 = I⃗[2]
    W1 = C1 / α * (I1 - 3)^α + C2 * (I1 - 3) + C3 / ζ * (I1 - 3)^ζ
    W2 = K1 / β * (I2 - 3)^β
    W = W1 + W2
    return W
end

function parameters(::Beda)
    return (:C1, :C2, :C3, :K1, :α, :β, :ζ)
end

function parameter_bounds(::Beda, data::AbstractHyperelasticTest)
    lb = (C1 = -Inf, C2 = -Inf, C3 = -Inf, K1 = -Inf, α = 0.0, β = 0.0, ζ = 1.0)
    ub = (C1 = Inf, C2 = Inf, C3 = Inf, K1 = Inf, α = 1.0, β = 1.0, ζ = Inf)
    return (lb = lb, ub = ub)
end

struct Amin{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = C_1 (I_1 - 3) + \\frac{C_2}{N + 1} (I_1 - 3)^{N + 1} + \\frac{C_3}{M + 1} (I_1 - 3)^{M + 1} + C_4 (I_2 - 3)
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `C1`
- `C2`
- `C3`
- `C4`
- `N`
- `M`

> Amin AF, Wiraguna SI, Bhuiyan AR, Okui Y. Hyperelasticity model for finite element analysis of natural and high damping rubbers in compression and shear. Journal of engineering mechanics. 2006 Jan;132(1):54-64.
"""
Amin(type::T = PrincipalValueForm()) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    Amin{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Amin{T},
    λ⃗::Vector{S},
    (; C1, C2, C3, C4, N, M),
) where {T<:PrincipalValueForm,S}
    return C1 * (I₁(λ⃗) - 3) +
           C2 / (N + 1) * (I₁(λ⃗) - 3)^(N + 1) +
           C3 / (M + 1) * (I₁(λ⃗) - 3)^(M + 1) +
           C4 * (I₂(λ⃗) - 3)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Amin{T},
    I⃗::Vector{S},
    (; C1, C2, C3, C4, N, M),
) where {T<:InvariantForm,S}
    return C1 * (I⃗[1] - 3) +
           C2 / (N + 1) * (I⃗[1] - 3)^(N + 1) +
           C3 / (M + 1) * (I⃗[1] - 3)^(M + 1) +
           C4 * (I⃗[2] - 3)
end

function parameters(::Amin)
    return (:C1, :C2, :C3, :C4, :N, :M)
end

struct LopezPamies{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = \\frac{3^{1 - \\alpha_i}}{2\\alpha_i} \\mu_i (I_1^{\\alpha_i} - 3^{\\alpha_i})
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `α⃗`
- `μ⃗`

> Lopez-Pamies O. A new I1-based hyperelastic model for rubber elastic materials. Comptes Rendus Mecanique. 2010 Jan 1;338(1):3-11.
"""
LopezPamies(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = LopezPamies{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::LopezPamies{T},
    λ⃗::Vector{S},
    (; α⃗, μ⃗),
) where {T<:PrincipalValueForm,S}
    @assert length(α⃗) == length(μ⃗) "length of α⃗ is not equal to length of μ⃗"
    I1 = I₁(λ⃗)
    return sum(@. (3^(1 - α⃗)) / (2α⃗) * μ⃗ * (I1^(α⃗) - 3^(α⃗)))
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::LopezPamies{T},
    I⃗::Vector{S},
    (; α⃗, μ⃗),
) where {T<:InvariantForm,S}
    @assert length(α⃗) == length(μ⃗) "length of α⃗ is not equal to length of μ⃗"
    return sum(@. (3^(1 - α⃗)) / (2α⃗) * μ⃗ * (I⃗[1]^(α⃗) - 3^(α⃗)))
end

function parameters(::LopezPamies)
    return (:α⃗, :μ⃗)
end

struct GenYeoh{T} <: AbstractIncompressibleModel{T} end

"""
$(SIGNATURES)

# Model:

```math
W = K_1 (I_1 - 3)^m + K_2 * (I_1 - 3)^p + K_3 * (I_1 - 3)^q
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `K1`
- `K2`
- `K3`
- `m`
- `p`
- `q`

> Hohenberger TW, Windslow RJ, Pugno NM, Busfield JJ. A constitutive model for both low and high strain nonlinearities in highly filled elastomers and implementation with user-defined material subroutines in ABAQUS. Rubber Chemistry and Technology. 2019;92(4):653-86.
"""
GenYeoh(type::T = PrincipalValueForm()) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    GenYeoh{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::GenYeoh{T},
    λ⃗::Vector{S},
    (; K1, K2, K3, m, p, q),
) where {T<:PrincipalValueForm,S}
    return K1 * (I₁(λ⃗) - 3)^m + K2 * (I₁(λ⃗) - 3)^p + K3 * (I₁(λ⃗) - 3)^q
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::GenYeoh{T},
    I⃗::Vector{S},
    (; K1, K2, K3, m, p, q),
) where {T<:InvariantForm,S}
    return K1 * (I⃗[1] - 3)^m + K2 * (I⃗[1] - 3)^p + K3 * (I⃗[1] - 3)^q
end

function parameters(::GenYeoh)
    return (:K1, :K2, :K3, :m, :p, :q)
end

struct HartSmith{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = \\frac{G\\exp{(-9k_1+k_1I_1)}}{k_1}+Gk_2\\log{I_2}
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `G`
- `k₁`
- `k₂`

> Hart-Smith LJ. Elasticity parameters for finite deformations of rubber-like materials. Zeitschrift für angewandte Mathematik und Physik ZAMP. 1966 Sep;17(5):608-26.
"""
HartSmith(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = HartSmith{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::HartSmith{T},
    λ⃗::Vector{S},
    (; G, k₁, k₂),
) where {T<:PrincipalValueForm,S}
    return G * exp(-9k₁ + k₁ * I₁(λ⃗)) / k₁ + G * k₂ * log(I₂(λ⃗))
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::HartSmith{T},
    I⃗::Vector{S},
    (; G, k₁, k₂),
) where {T<:InvariantForm,S}
    return G * exp(-9k₁ + k₁ * I⃗[1]) / k₁ + G * k₂ * log(I⃗[2])
end

function parameters(::HartSmith)
    return (:G, :k₁, :k₂)
end

struct VerondaWestmann{T} <: AbstractIncompressibleModel{T} end

"""
$(SIGNATURES)

# Model:

```math
W = C_1 (\\exp(\\alpha(I_1 - 3)) - 1) + C_2 (I_2 - 3)
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `C1`
- `C2`
- `α`

> Veronda DR, Westmann RA. Mechanical characterization of skin—finite deformations. Journal of biomechanics. 1970 Jan 1;3(1):111-24.
"""
VerondaWestmann(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = VerondaWestmann{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::VerondaWestmann{T},
    λ⃗::Vector{S},
    (; C1, C2, α),
) where {T<:PrincipalValueForm,S}
    return C1 * (exp(α * (I₁(λ⃗) - 3)) - 1) + C2 * (I₂(λ⃗) - 3)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::VerondaWestmann{T},
    I⃗::Vector{S},
    (; C1, C2, α),
) where {T<:InvariantForm,S}
    return C1 * (exp(α * (I⃗[1] - 3)) - 1) + C2 * (I⃗[2] - 3)
end

function parameters(::VerondaWestmann)
    return (:C1, :C2, :α)
end

struct FungDemiray{T} <: AbstractIncompressibleModel{T} end

"""
$(SIGNATURES)

# Model:

```math
W = \\frac{\\mu}{2 * b} (\\exp(b(I_1 - 3)) - 1)
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `μ`
- `b`

> Fung YC. Elasticity of soft tissues in simple elongation. American Journal of Physiology-Legacy Content. 1967 Dec 1;213(6):1532-44.
> Demiray H. A note on the elasticity of soft biological tissues. Journal of biomechanics. 1972 May 1;5(3):309-11.
"""
FungDemiray(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = FungDemiray{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::FungDemiray{T},
    λ⃗::Vector{S},
    (; μ, b),
) where {T<:PrincipalValueForm,S}
    return μ / (2 * b) * (exp(b * (I₁(λ⃗) - 3)) - 1)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::FungDemiray{T},
    I⃗::Vector{S},
    (; μ, b),
) where {T<:InvariantForm,S}
    return μ / (2 * b) * (exp(b * (I⃗[1] - 3)) - 1)
end

function parameters(::FungDemiray)
    return (:μ, :b)
end

struct Vito{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = \\alpha (\\exp\\bigg(\\beta (I_1 - 3)\\bigg) + \\gamma  (I_2 - 3)) - 1)
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `α`
- `β`
- `γ`

> Vito R. A note on arterial elasticity. Journal of Biomechanics. 1973 Sep 1;6(5):561-4.
"""
Vito(type::T = PrincipalValueForm()) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    Vito{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Vito{T},
    λ⃗::Vector{S},
    (; α, β, γ),
) where {T<:PrincipalValueForm,S}
    return α * (exp(β * (I₁(λ⃗) - 3) + γ * (I₂(λ⃗) - 3)) - 1)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Vito{T},
    I⃗::Vector{S},
    (; α, β, γ),
) where {T<:InvariantForm,S}
    return α * (exp(β * (I⃗[1] - 3) + γ * (I⃗[2] - 3)) - 1)
end

function parameters(::Vito)
    return (:α, :β, :γ)
end

struct ModifiedYeoh{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = C_{10} * (I_1 - 3) + C_{20} * (I_1 - 3)^2 + C_{30} * (I_1 - 3)^3 + \\alpha / \\beta * (1 - \\exp{-\\beta * (I_1 - 3)})
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `C10`
- `C20`
- `C30`
- `α`
- `β`

> He H, Zhang Q, Zhang Y, Chen J, Zhang L, Li F. A comparative study of 85 hyperelastic constitutive models for both unfilled rubber and highly filled rubber nanocomposite material. Nano Materials Science. 2021 Jul 16.
"""
ModifiedYeoh(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = ModifiedYeoh{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::ModifiedYeoh{T},
    λ⃗::Vector{S},
    (; C10, C20, C30, α, β),
) where {T<:PrincipalValueForm,S}
    return C10 * (I₁(λ⃗) - 3) +
           C20 * (I₁(λ⃗) - 3)^2 +
           C30 * (I₁(λ⃗) - 3)^3 +
           α / β * (1 - exp(-β * (I₁(λ⃗) - 3)))
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::ModifiedYeoh{T},
    I⃗::Vector{S},
    (; C10, C20, C30, α, β),
) where {T<:InvariantForm,S}
    return C10 * (I⃗[1] - 3) +
           C20 * (I⃗[1] - 3)^2 +
           C30 * (I⃗[1] - 3)^3 +
           α / β * (1 - exp(-β * (I⃗[1] - 3)))
end

function parameters(::ModifiedYeoh)
    return (:C10, :C20, :C30, :α, :β)
end

struct ChevalierMarco{T} <: AbstractIncompressibleModel{T}
    ∂W∂I1::Function
    ∂W∂I2::Function
end

"""
$(SIGNATURES)

# Model:

```math
W = \\int\\limits_{3}^{I_1(\\vec\\lambda)} \\exp\\bigg(\\sum\\limits_{i=0}^{N}a_i(I_1-3)^i\\bigg)\\text{d}I_1+ \\int\\limits_{3}^{I_2(\\vec\\lambda)} \\sum\\limits_{i=0}^{n}\\frac{b_i}{I_2^i}\\text{d}I_2
```

```math
[\\mathbf{S}] = 2(I-\\frac{\\partial W}{\\partial I_1} - C^{-2}\\frac{\\partial W}{\\partial I_2})
```

```math
[\\mathbf{\\sigma}] = \\mathbf{F} \\cdot \\mathbf{S}
```

# Parameters:
- `a⃗`
- `b⃗`

Note:
- Model is not compatible with AD. A method for accessing the Second Piola Kirchoff Tensor and Cauchy Stress Tensor have been implemented.

> Chevalier L, Marco Y. Tools for multiaxial validation of behavior laws chosen for modeling hyper‐elasticity of rubber‐like materials. Polymer Engineering & Science. 2002 Feb;42(2):280-98.
"""
function ChevalierMarco(::T = PrincipalValueForm()) where {T<:Union{PrincipalValueForm}}
    function ∂W∂I1(I₁, a⃗)
        L_a = size(a⃗, 1)
        return exp(sum(@. a⃗ * (I₁ - 3)^(1:L_a)))
    end
    function ∂W∂I2(I₂, b⃗)
        L_b = size(b⃗, 1)
        return sum(@. b⃗ / I₂^(1:L_b))
    end
    ChevalierMarco{T}(∂W∂I1, ∂W∂I2)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    W::ChevalierMarco{T},
    λ⃗::Vector{S},
    (; a⃗, b⃗),
) where {T<:PrincipalValueForm,S}
    return quadgk(Base.Fix2(W.∂W∂I1, a⃗), 3, I₁(λ⃗))[1] +
           quadgk(Base.Fix2(W.∂W∂I2, b⃗), 3, I₂(λ⃗))[1]
end

# function ContinuumMechanicsBase.StrainEnergyDensity(W::ChevalierMarco{T}, I⃗::Vector{S}, (; a⃗, b⃗)) where {T<:InvariantForm, S}
#     # ∂W∂I1(I₁) = exp(sum(@tullio _ := a⃗[i] * (I₁ - 3)^(i - 1)))
#     # ∂W∂I2(I₂) = @tullio _ := b⃗[i] / I₂^(i - 1)
#     return quadgk(Base.Fix2(W.∂W∂I1,a⃗), 3, I⃗[1])[1] + quadgk(Base.Fix2(W.∂W∂I2,b⃗), 3, I⃗[2])[1]
# end

function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    W::ChevalierMarco{T},
    λ⃗::Vector{S},
    (; a⃗, b⃗);
    kwargs...,
) where {T<:PrincipalValueForm,S}
    𝐒 = 2 * (I(3) * W.∂W∂I1(I₁(λ⃗), a⃗) - diagm(λ⃗ .^ 2)^(-2) * W.∂W∂I2(I₂(λ⃗), b⃗))
    sᵢ = diag(𝐒)
    sᵢ = sᵢ
    return sᵢ
end

function ContinuumMechanicsBase.CauchyStressTensor(
    W::ChevalierMarco{T},
    λ⃗::Vector{S},
    p;
    kwargs...,
) where {T<:PrincipalValueForm,S}
    s = ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(W, λ⃗, p)
    σ = λ⃗ .* s
    return σ
end

function parameters(::ChevalierMarco)
    return (:a⃗, :b⃗)
end

struct GornetDesmorat{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:
```math
W = h_1\\int\\exp{h_3(I_1-3)^2}\\text{d}I_1+3h_2\\int\\frac{1}{\\sqrt{I_2}}\\text{d}I_2 = \\frac{h_1 \\sqrt{\\pi} \\text{erfi}(\\sqrt{h_3}(I_1-3)^2)}{2\\sqrt{h_3}}+6h_2\\sqrt{I_2}
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- h₁
- h₂
- h₃

# Notes:
- The differential form was original form and the closed form SEF was determine via symbolic integration in Mathematica. The model is not compatible with AD and has methods for the Second Piola Kirchoff Stress Tensor and Cauchy Stress Tensor implemented.

> Gornet L, Marckmann G, Desmorat R, Charrier P. A new isotropic hyperelastic strain energy function in terms of invariants and its derivation into a pseudo-elastic model for Mullins effect: application to finite element analysis. Constitutive Models for Rubbers VII. 2012:265-71.
"""
GornetDesmorat(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = GornetDesmorat{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::GornetDesmorat{T},
    λ⃗::Vector{S},
    (; h₁, h₂, h₃),
) where {T<:PrincipalValueForm,S}
    return h₁ * √π * erfi(√h₃ * (I₁(λ⃗) - 3)^2) / 2 / √h₃ + 6 * h₂ * √(I₂(λ⃗))
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::GornetDesmorat{T},
    I⃗::Vector{S},
    (; h₁, h₂, h₃),
) where {T<:InvariantForm,S}
    return h₁ * √π * erfi(√h₃ * (I⃗[1] - 3)^2) / 2 / √h₃ + 6 * h₂ * √(I⃗[2])
end

function ContinuumMechanicsBase.CauchyStressTensor(
    ::GornetDesmorat{T},
    λ⃗::Vector{S},
    (; h₁, h₂, h₃);
    kwargs...,
) where {T<:PrincipalValueForm,S}
    B = λ⃗ .^ 2
    _I₁ = I₁(λ⃗)
    _I₂ = I₂(λ⃗)
    ∂W∂I₁ = h₁ * exp(h₃ * (_I₁ - 3)^2)
    ∂W∂I₂ = 3 * h₂ * exp(1 / sqrt(_I₂))
    σ = 2 * (∂W∂I₁ + _I₁ * ∂W∂I₂) * B - 2 * ∂W∂I₂ * (B .^ 2)
    return σ
end

function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    ψ::GornetDesmorat{T},
    λ⃗::Vector{S},
    ps;
    kwargs...,
) where {T<:PrincipalValueForm,S}
    σ = CauchyStressTensor(ψ, λ⃗, ps; kwargs...)
    s = σ ./ λ⃗
    return s
end

function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    ψ::GornetDesmorat{T},
    F::Matrix{R},
    (; h₁, h₂, h₃);
    ad_type = nothing,
    kwargs...,
) where {T<:InvariantForm,R}
    I1 = I₁(F)
    I2 = I₂(F)
    I3 = I₃(F)
    ∂W∂I₁ = h₁ * exp(h₃ * (I1 - 3)^2)
    ∂W∂I₂ = 3 * h₂ * exp(1 / sqrt(I2))
    ∂ψ∂I = [∂W∂I₁, ∂W∂I₂, 0.0]
    S = 2∂ψ∂I[1] * F' + 2∂ψ∂I[2] * (I1 * F' + F' * F * F') + 2I3 * ∂ψ∂I[3] * inv(F)
    return S
end

function ContinuumMechanicsBase.CauchyStressTensor(
    ψ::GornetDesmorat{T},
    F::Matrix{S},
    (; h₁, h₂, h₃);
    ad_type,
    kwargs...,
) where {T<:InvariantForm,S}
    I1 = I₁(F)
    I2 = I₂(F)
    I3 = I₃(F)
    J = sqrt(I3)
    ∂W∂I₁ = h₁ * exp(h₃ * (I1 - 3)^2)
    ∂W∂I₂ = 3 * h₂ * exp(1 / sqrt(I2))
    ∂ψ∂I = [∂W∂I₁, ∂W∂I₂, 0.0]
    B = F * F'
    σ =
        2 * inv(J) * (∂ψ∂I[1] + I1 * ∂ψ∂I[2]) * B - 2 * inv(J) * ∂ψ∂I[2] * B^2 +
        2 * J * ∂ψ∂I[3] * I
    return σ
end


function parameters(::GornetDesmorat)
    return (:h₁, :h₂, :h₃)
end

struct MansouriDarijani{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = A_1\\exp{m_1(I_1-3)-1}+B_1\\exp{n_1(I_2-3)-1}
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `A1`
- `m1`
- `B1`
- `n1`

> Mansouri MR, Darijani H. Constitutive modeling of isotropic hyperelastic materials in an exponential framework using a self-contained approach. International Journal of Solids and Structures. 2014 Dec 1;51(25-26):4316-26.
"""
MansouriDarijani(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = MansouriDarijani{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::MansouriDarijani{T},
    λ⃗::Vector{S},
    (; A1, m1, B1, n1),
) where {T<:PrincipalValueForm,S}
    return A1 * (exp(m1 * (I₁(λ⃗) - 3)) - 1) + B1 * (exp(n1 * (I₂(λ⃗) - 3)) - 1)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::MansouriDarijani{T},
    I⃗::Vector{S},
    (; A1, m1, B1, n1),
) where {T<:InvariantForm,S}
    return A1 * (exp(m1 * (I⃗[1] - 3)) - 1) + B1 * (exp(n1 * (I⃗[2] - 3)) - 1)
end

function parameters(::MansouriDarijani)
    return (:A1, :m1, :B1, :n1)
end

struct GentThomas{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = C_1(I_1-3)+C_2\\log(\\frac{I_2}{3})
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Paramters:
- `C1`
- `C2`

> Gent AN, Thomas AG. Forms for the stored (strain) energy function for vulcanized rubber. Journal of Polymer Science. 1958 Apr;28(118):625-8.
"""
GentThomas(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = GentThomas{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::GentThomas{T},
    λ⃗::Vector{S},
    (; C1, C2),
) where {T<:PrincipalValueForm,S}
    return C1 * (I₁(λ⃗) - 3) + C2 * log(I₂(λ⃗) / 3)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::GentThomas{T},
    I⃗::Vector{S},
    (; C1, C2),
) where {T<:InvariantForm,S}
    return C1 * (I⃗[1] - 3) + C2 * log(I⃗[2] / 3)
end

function parameters(::GentThomas)
    return (:C1, :C2)
end

struct LambertDianiRey{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = \\int\\limits_{3}^{I_1}\\exp\\bigg(\\sum\\limits_{i=0}^{n}a_i(I_1-3)^i\\bigg)\\text{d}I_1+\\int\\limits_{3}^{I_2}\\sum\\limits_{j=0}^{m}b_i\\log(I_2)^i\\text{d}I_2
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `a⃗`
- `b⃗`

> Lambert-Diani J, Rey C. New phenomenological behavior laws for rubbers and thermoplastic elastomers. European Journal of Mechanics-A/Solids. 1999 Nov 1;18(6):1027-43.
"""
LambertDianiRey(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = LambertDianiRey{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::LambertDianiRey{T},
    λ⃗::Vector{S},
    (; a⃗, b⃗),
) where {T<:PrincipalValueForm,S}
    length_a = length(a⃗)
    length_b = length(b⃗)
    ∂W∂I₁(I1) = exp(sum(@. a⃗ * (I1 - 3)^(1:length_a)))
    ∂W∂I₂(I2) = exp(sum(@. b⃗ * (I2 - 3)^(1:length_b)))
    # @. b⃗*(I2-3)^(1:length_b)
    # ∂W∂I₁(I₁) = exp(@tullio _ := a⃗[i] .* (I₁ .- 3) .^ i)
    # ∂W∂I₂(I₂) = exp(@tullio _ := b⃗[i] .* log(I₂) .^ i)
    return quadgk(∂W∂I₁, 3, I₁(λ⃗))[1] + quadgk(∂W∂I₂, 3, I₂(λ⃗))[1]
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::LambertDianiRey{T},
    I⃗::Vector{S},
    (; a⃗, b⃗),
) where {T<:InvariantForm,S}
    # ∂W∂I₁(I₁) = exp(@tullio _ := a⃗[i] .* (I₁ .- 3) .^ i)
    # ∂W∂I₂(I₂) = exp(@tullio _ := b⃗[i] .* log(I₂) .^ i)
    length_a = length(a⃗)
    length_b = length(b⃗)
    ∂W∂I₁(I1) = exp(sum(@. a⃗ * (I1 - 3)^(1:length_a)))
    ∂W∂I₂(I2) = exp(sum(@. b⃗ * (I2 - 3)^(1:length_b)))
    return quadgk(∂W∂I₁, 3, I⃗[1])[1] + quadgk(∂W∂I₂, 3, I⃗[2])[1]
end


function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    ::LambertDianiRey{T},
    λ⃗::Vector{S},
    (; a⃗, b⃗);
    kwargs...,
) where {T<:PrincipalValueForm,S}
    # ∂W∂I₁ = exp(@tullio _ := a⃗[i] .* (I₁(λ⃗) .- 3) .^ i)
    # ∂W∂I₂ = exp(@tullio _ := b⃗[i] .* log(I₂(λ⃗)) .^ i)
    length_a = length(a⃗)
    length_b = length(b⃗)
    I1 = I₁(λ⃗)
    I2 = I₂(λ⃗)
    ∂W∂I₁ = exp(sum(@. a⃗ * (I1 - 3)^(1:length_a)))
    ∂W∂I₂ = exp(sum(@. b⃗ * (I2 - 3)^(1:length_b)))
    𝐒 = 2 * (I * ∂W∂I₁ - diagm(λ⃗ .^ 2)^(-2) * ∂W∂I₂)
    sᵢ = diag(𝐒)
    # sᵢ = sᵢ .- sᵢ[3] .* λ⃗[3] ./ λ⃗
    return sᵢ
end

function ContinuumMechanicsBase.CauchyStressTensor(
    ψ::LambertDianiRey{T},
    λ⃗::Vector{S},
    ps;
    kwargs...,
) where {T<:PrincipalValueForm,S}
    s = ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(ψ, λ⃗, ps)
    σᵢ = λ⃗ .* s
    return σᵢ
end

function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    ::LambertDianiRey{T},
    F::Matrix{R},
    p;
    kwargs...,
) where {T<:InvariantForm,R}
    (; a⃗, b⃗) = p
    I⃗ = [I₁(F), I₂(F), I₃(F)]
    length_a = length(a⃗)
    length_b = length(b⃗)
    ∂W∂I₁ = exp(sum(@. a⃗ * (I⃗[1] - 3)^(1:length_a)))
    ∂W∂I₂ = exp(sum(@. b⃗ * (I⃗[2] - 3)^(1:length_b)))
    ∂W∂I₃ = zero(eltype(I⃗))
    ∂ψ∂I = [∂W∂I₁, ∂W∂I₂, ∂W∂I₃]
    S = 2∂ψ∂I[1] * F' + 2∂ψ∂I[2] * (I⃗[1] * F' + F' * F * F') + 2 * I⃗[3] * ∂ψ∂I[3] * inv(F)
    return S
end

function ContinuumMechanicsBase.CauchyStressTensor(
    ψ::LambertDianiRey{T},
    F::Matrix{S},
    p;
    kwargs...,
) where {T<:InvariantForm,S}
    (; a⃗, b⃗) = p
    I⃗ = [I₁(F), I₂(F), I₃(F)]
    J = sqrt(I₃(F))
    B = F * F'
    (; a⃗, b⃗) = p
    length_a = length(a⃗)
    length_b = length(b⃗)
    ∂W∂I₁ = exp(sum(@. a⃗ * (I⃗[1] - 3)^(1:length_a)))
    ∂W∂I₂ = exp(sum(@. b⃗ * (I⃗[2] - 3)^(1:length_b)))
    ∂W∂I₃ = zero(eltype(I⃗))
    ∂ψ∂I = [∂W∂I₁, ∂W∂I₂, ∂W∂I₃]
    σ =
        2 * inv(J) * (∂ψ∂I[1] + I⃗[1] * ∂ψ∂I[2]) * B - 2 * inv(J) * ∂ψ∂I[2] * B^2 +
        2 * sqrt(I⃗[3]) * ∂ψ∂I[3] * I
    return σ
end


function parameters(::LambertDianiRey)
    return (:a⃗, :b⃗)
end

struct HossMarczakI{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = \\frac{\\alpha}{\\beta}(1-\\exp{-\\beta(I_1-3)})+\\frac{\\mu}{2b}\\bigg((1+\\frac{b}{n}(I_1-3))^n -1\\bigg)
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `α`
- `β`
- `μ`
- `b`
- `n`

# Note:
- The authors suggested this model for low strains.

> Hoss L, Marczak RJ. A new constitutive model for rubber-like materials. Mecánica Computacional. 2010;29(28):2759-73.
"""
HossMarczakI(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = HossMarczakI{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::HossMarczakI{T},
    λ⃗::Vector{S},
    (; α, β, μ, b, n),
) where {T<:PrincipalValueForm,S}
    return α / β * (1 - exp(-β * (I₁(λ⃗) - 3))) +
           μ / (2b) * ((1 + b / n * (I₁(λ⃗) - 3))^n - 1)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::HossMarczakI{T},
    I⃗::Vector{S},
    (; α, β, μ, b, n),
) where {T<:InvariantForm,S}
    return α / β * (1 - exp(-β * (I⃗[1] - 3))) +
           μ / (2b) * ((1 + b / n * (I⃗[1] - 3))^n - 1)
end

function parameters(::HossMarczakI)
    return (:α, :β, :μ, :b, :n)
end

function parameter_bounds(::HossMarczakI, data::AbstractHyperelasticTest)
    lb = (α = -Inf, β = 0, μ = -Inf, b = 0, n = 0)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct HossMarczakII{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = \\frac{\\alpha}{\\beta}(1-\\exp{-\\beta(I_1-3)})+\\frac{\\mu}{2b}\\bigg((1+\\frac{b}{n}(I_1-3))^n -1\\bigg)+C_2\\log(\\frac{I_2}{3})
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `α`
- `β`
- `μ`
- `b`
- `n`
- `C2`

# Note:
- The authors suggests this model for high strains.

> Hoss L, Marczak RJ. A new constitutive model for rubber-like materials. Mecánica Computacional. 2010;29(28):2759-73.
"""
HossMarczakII(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = HossMarczakII{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::HossMarczakII{T},
    λ⃗::Vector{S},
    (; α, β, μ, b, n, C2),
) where {T<:PrincipalValueForm,S}
    return α / β * (1 - exp(-β * (I₁(λ⃗) - 3))) +
           μ / (2b) * ((1 + b / n * (I₁(λ⃗) - 3))^n - 1) +
           C2 * log(I₂(λ⃗) / 3)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::HossMarczakII{T},
    I⃗::Vector{S},
    (; α, β, μ, b, n, C2),
) where {T<:InvariantForm,S}
    return α / β * (1 - exp(-β * (I⃗[1] - 3))) +
           μ / (2b) * ((1 + b / n * (I⃗[1] - 3))^n - 1) +
           C2 * log(I⃗[2] / 3)
end

function parameters(::HossMarczakII)
    return (:α, :β, :μ, :b, :n, :C2)
end

function parameter_bounds(::HossMarczakII, data::AbstractHyperelasticTest)
    lb = (α = -Inf, β = 0, μ = -Inf, b = 0, n = 0, C2 = -Inf)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct ExpLn{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = A\\bigg[\\frac{1}{a}\\exp{(a(I_1-3))}+b(I_1-2)(1-\\log{I_1-2})-\\frac{1}{a}-b\\bigg]
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `A`
- `a`
- `b`

> Khajehsaeid H, Arghavani J, Naghdabadi R. A hyperelastic constitutive model for rubber-like materials. European Journal of Mechanics-A/Solids. 2013 Mar 1;38:144-51.
"""
ExpLn(type::T = PrincipalValueForm()) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    ExpLn{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::ExpLn{T},
    λ⃗::Vector{S},
    (; A, a, b),
) where {T<:PrincipalValueForm,S}
    return A * (
        1 / a * exp(a * (I₁(λ⃗) - 3)) + b * (I₁(λ⃗) - 2) * (1 - log(I₁(λ⃗) - 2)) - 1 / a -
        b
    )
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::ExpLn{T},
    I⃗::Vector{S},
    (; A, a, b),
) where {T<:InvariantForm,S}
    return A * (
        1 / a * exp(a * (I⃗[1] - 3)) + b * (I⃗[1] - 2) * (1 - log(I⃗[1] - 2)) - 1 / a - b
    )
end

function parameters(::ExpLn)
    return (:A, :a, :b)
end

struct VanDerWaals{T} <: AbstractIncompressibleModel{T} end

"""
$(SIGNATURES)

# Model:

```math
W = -\\mu\\{(\\lambda_m^2-3)\\log(1-\\Theta)+\\Theta\\}-\\frac{2\\alpha}{3}\\bigg(\\frac{I-3}{2}\\bigg)^{3/2}
```

where:

```math
\\Theta = \\frac{\\beta I_1 + (1-\\beta)I_2-3}{\\lambda_m^2-3)}
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- μ
- λm
- β
- α

> Kilian HG, Enderle HF, Unseld K. The use of the van der Waals model to elucidate universal aspects of structure-property relationships in simply extended dry and swollen rubbers. Colloid and Polymer Science. 1986 Oct;264(10):866-76.
> Ambacher H, Enderle HF, Kilian HG, Sauter A. Relaxation in permanent networks. InRelaxation in Polymers 1989 (pp. 209-220). Steinkopff.
> Kilian HG. A molecular interpretation of the parameters of the van der Waals equation of state for real networks. Polymer Bulletin. 1980 Sep;3(3):151-8.
"""
VanDerWaals(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = VanDerWaals{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::VanDerWaals{T},
    λ⃗::Vector{S},
    (; μ, λm, β, α),
) where {T<:PrincipalValueForm,S}
    I = β * I₁(λ⃗) + (1 - β) * I₂(λ⃗)
    θ = sqrt((I - 3) / (λm^2 - 3))
    W = -μ * ((λm^2 - 3) * log(1 - θ) + θ) - (2 * α / 3) * sqrt((I - 3) / 2)^3
    return W
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::VanDerWaals{T},
    I⃗::Vector{S},
    (; μ, λm, β, α),
) where {T<:InvariantForm,S}
    I = β * I⃗[1] + (1 - β) * I⃗[2]
    θ = (I - 3) / (λm^2 - 3)
    return μ * (-(λm^2 - 3) * log(1 - θ) + θ) - 2 / 3 * α * ((I - 3) / 2)^(3 / 2)
end

function parameter_bounds(::VanDerWaals, data::AbstractHyperelasticTest)
    _I2 = data.data.λ .|> I₂
    _I1 = data.data.λ .|> I₁
    β_min = maximum(@. (3 - _I2) / (_I1 - _I2))
    lb = (μ = 0.0, λm = sqrt(3), β = β_min, α = 0.0)
    ub = (μ = Inf, λm = Inf, β = 1.0, α = Inf)
    return (ub = ub, lb = lb)
end

function parameters(::VanDerWaals)
    return (:μ, :λm, :β, :α)
end

# function constraints(::VanDerWaals, data::AbstractHyperelasticTest)
#     I₁_max = maximum(I₁.(data.data.λ))
#     I₂_max = maximum(I₂.(data.data.λ))
#     return f(u, p) = [1 - (u.β * I₁_max + (1 - u.β) * I₂_max - 3) / (u.λm^2 - 3)]
# end

struct Gent{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = -\\frac{\\mu J_m}{2}\\log{\\bigg(1-\\frac{I_1-3}{J_m}\\bigg)}
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `μ`:  Small strain shear modulus
- `Jₘ`: Limiting stretch invariant

> Gent AN. A new constitutive relation for rubber. Rubber chemistry and technology. 1996 Mar;69(1):59-61.
"""
Gent(type::T = PrincipalValueForm()) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    Gent{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Gent{T},
    λ⃗::Vector{S},
    p,
) where {T<:PrincipalValueForm,S}
    (; μ, Jₘ) = p
    return -(μ * Jₘ) / 2 * log(1 - (I₁(λ⃗) - 3) / Jₘ)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Gent{T},
    I⃗::Vector{S},
    (; μ, Jₘ),
) where {T<:InvariantForm,S}
    -(μ * Jₘ) / 2 * log(1 - (I⃗[1] - 3) / Jₘ)
end

function parameters(::Gent)
    return (:μ, :Jₘ)
end

function parameter_bounds(::Gent, test::AbstractHyperelasticTest{S,T}) where {S,T}
    I₁_max = maximum(I₁.(test.data.λ))
    Jₘ_min = I₁_max - 3
    lb = (μ = zero(T), Jₘ = Jₘ_min)
    ub = nothing
    return (lb = lb, ub = ub)
end


struct TakamizawaHayashi{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = -c\\log\\left1-\\left(\\frac{I_1-3}{J_m}\\right)^2\\right
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `c`
- `Jₘ`

> Takamizawa K, Hayashi K. Strain energy density function and uniform strain hypothesis for arterial mechanics. Journal of biomechanics. 1987 Jan 1;20(1):7-17.
"""
TakamizawaHayashi(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = TakamizawaHayashi{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::TakamizawaHayashi{T},
    λ⃗::Vector{S},
    (; c, Jₘ),
) where {T<:PrincipalValueForm,S}
    return -c * log(1 - ((I₁(λ⃗) - 3) / Jₘ)^2)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::TakamizawaHayashi{T},
    I⃗::Vector{S},
    (; c, Jₘ),
) where {T<:InvariantForm,S}
    return -c * log(1 - ((I⃗[1] - 3) / Jₘ)^2)
end

function parameters(::TakamizawaHayashi)
    return (:c, :Jₘ)
end

function parameter_bounds(::TakamizawaHayashi, data::AbstractHyperelasticTest)
    I₁_max = maximum(I₁.(data.data.λ))
    Jₘ_min = I₁_max - 3
    lb = (c = -Inf, Jₘ = Jₘ_min)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct YeohFleming{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = \\frac{A}{B}(1-\\exp{-B(I_1-3)}) - C_{10}(I_m-3)\\log{1-\\frac{I_1-3}{I_m-3}}
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `A`
- `B`
- `C10`
- `Im`

>  Yeoh OH, Fleming PD. A new attempt to reconcile the statistical and phenomenological theories of rubber elasticity. Journal of Polymer Science Part B: Polymer Physics. 1997 Sep 15;35(12):1919-31.
"""
YeohFleming(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = YeohFleming{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::YeohFleming{T},
    λ⃗::Vector{S},
    (; A, B, C10, Im),
) where {T<:PrincipalValueForm,S}
    return A / B * (1 - exp(-B * (I₁(λ⃗) - 3))) -
           C10 * (Im - 3) * log(1 - ((I₁(λ⃗) - 3) / (Im - 3)))
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::YeohFleming{T},
    I⃗::Vector{S},
    (; A, B, C10, Im),
) where {T<:InvariantForm,S}
    return A / B * (1 - exp(-B * (I⃗[1] - 3))) -
           C10 * (Im - 3) * log(1 - ((I⃗[1] - 3) / (Im - 3)))
end

function parameters(::YeohFleming)
    return (:A, :B, :C10, :Im)
end

function parameter_bounds(::YeohFleming, data::AbstractHyperelasticTest)
    Iₘ_min = maximum(I₁, data.data.λ)
    lb = (A = -Inf, B = -Inf, C10 = -Inf, Im = Iₘ_min)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct PucciSaccomandi{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = K\\log{\\frac{I_2}{3}}-\\frac{\\mu J_m}{2}\\log{1-\\frac{I_1-3}{J-m}}
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `K`
- `μ`
- `Jₘ`

> Pucci E, Saccomandi G. A note on the Gent model for rubber-like materials. Rubber chemistry and technology. 2002 Nov;75(5):839-52.
"""
PucciSaccomandi(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = PucciSaccomandi{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::PucciSaccomandi{T},
    λ⃗::Vector{S},
    (; K, μ, Jₘ),
) where {T<:PrincipalValueForm,S}
    return K * log(I₂(λ⃗) / 3) - μ * Jₘ / 2 * log(1 - (I₁(λ⃗) - 3) / Jₘ)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::PucciSaccomandi{T},
    I⃗::Vector{S},
    (; K, μ, Jₘ),
) where {T<:InvariantForm,S}
    return K * log(I⃗[2] / 3) - μ * Jₘ / 2 * log(1 - (I⃗[1] - 3) / Jₘ)
end

function parameters(::PucciSaccomandi)
    return (:K, :μ, :Jₘ)
end

function parameter_bounds(::PucciSaccomandi, data::AbstractHyperelasticTest)
    I₁_max = maximum(I₁.(data.data.λ))
    Jₘ_min = I₁_max - 3
    lb = (K = -Inf, μ = -Inf, Jₘ = Jₘ_min)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct HorganSaccomandi{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = -\\frac{\\mu J}{2}\\log\\bigg(\\frac{J^3-J^2I_1+JI_2-1}{(J-1)^3}\\bigg)
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `μ`
- `J`

> Horgan CO, Saccomandi G. Constitutive models for compressible nonlinearly elastic materials with limiting chain extensibility. Journal of Elasticity. 2004 Nov;77(2):123-38.\
> Horgan CO, Saccomandi G. Constitutive models for atactic elastomers. InWaves And Stability In Continuous Media 2004 (pp. 281-294).
"""
HorganSaccomandi(
    type::T = PrincipalValueForm(),
) where {T<:Union{InvariantForm,PrincipalValueForm}} = HorganSaccomandi{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::HorganSaccomandi{T},
    λ⃗::Vector{S},
    (; μ, J),
) where {T<:PrincipalValueForm,S}
    return -μ * J / 2 * log((J^3 - J^2 * I₁(λ⃗) + J * I₂(λ⃗) - 1) / (J - 1)^3)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::HorganSaccomandi{T},
    I⃗::Vector{S},
    (; μ, J),
) where {T<:InvariantForm,S}
    return -μ * J / 2 * log((J^3 - J^2 * I⃗[1] + J * I⃗[2] - 1) / (J - 1)^3)
end

function parameters(::HorganSaccomandi)
    return (:μ, :J)
end

function parameter_bounds(::HorganSaccomandi, data::AbstractHyperelasticTest)
    _I1 = @. I₁(data.data.λ)
    _I2 = @. I₂(data.data.λ)

    Js = @. 1 / 6 * (
        2 * _I1 +
        (2 * (2^(1 / 3)) * (_I1^2 - 3 * _I2)) / cbrt(27 + 2 * (_I1^3) - 9 * _I1 * _I2) +
        2^(2 / 3) * cbrt(27 + 2 * (_I1^3) - 9 * _I1 * _I2)
    )

    J_min = maximum(Js[(!isnan).(Js)])

    lb = (μ = -Inf, J = J_min)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct Beatty{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = -\\frac{G_0 I_m(I_m-3)}{2(2I_m-3)}\\log\\bigg(\\frac{1-\\frac{I_1-3}{I_m-3}}{1+\\frac{I_1-3}{I_m}} \\bigg)
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- `G₀`
- `Iₘ`

> Beatty MF. On constitutive models for limited elastic, molecular based materials. Mathematics and mechanics of solids. 2008 Jul;13(5):375-87.
"""
Beatty(type::T = PrincipalValueForm()) where {T<:Union{InvariantForm,PrincipalValueForm}} =
    Beatty{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Beatty{T},
    λ⃗::Vector{S},
    (; G₀, Iₘ),
) where {T<:PrincipalValueForm,S}
    return -G₀ * Iₘ * (Iₘ - 3) / 2 / (2Iₘ - 3) *
           log((1 - (I₁(λ⃗) - 3) / (Iₘ - 3)) / (1 + (I₁(λ⃗) - 3) / (Iₘ)))
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Beatty{T},
    I⃗::Vector{S},
    (; G₀, Iₘ),
) where {T<:InvariantForm,S}
    return -G₀ * Iₘ * (Iₘ - 3) / 2 / (2Iₘ - 3) *
           log((1 - (I⃗[1] - 3) / (Iₘ - 3)) / (1 + (I⃗[1] - 3) / (Iₘ)))
end

function parameters(::Beatty)
    return (:G₀, :Iₘ)
end

function parameter_bounds(::Beatty, data::AbstractHyperelasticTest)
    Iₘ_min = maximum(I₁, data.data.λ)
    lb = (G₀ = -Inf, Iₘ = Iₘ_min)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct HorganMurphy{T} <: AbstractIncompressibleModel{T} end

"""
$(SIGNATURES)

# Model:

```math
W = -\\frac{2\\mu J_m}{c^2}\\log\\left(1-\\frac{\\lambda_1^c+\\lambda_2^c+\\lambda_3^c-3}{J_m}\\right)
```

# Parameters:
- `μ`
- `Jₘ`
- `c`

> Horgan CO, Murphy JG. Limiting chain extensibility constitutive models of Valanis–Landel type. Journal of Elasticity. 2007 Feb;86(2):101-11.
"""
HorganMurphy() = HorganMurphy{PrincipalValueForm}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::HorganMurphy{T},
    λ⃗::Vector{S},
    (; μ, Jₘ, c),
) where {T<:PrincipalValueForm,S}
    return -2 * μ * Jₘ / c^2 * log(1 - (sum(λ⃗ .^ c) - 3) / Jₘ)
    # -2 * ps.μ  * ps.J / ps.c^2 * log(1 - (sum(λ⃗ .^ ps.c) - 3) / ps.J)
end

function parameters(::HorganMurphy)
    return (:μ, :Jₘ, :c)
end

# function constraints(::HorganMurphy, data::AbstractHyperelasticTest)
#     function f(res, u, p)
#         max_sum = minimum(λ⃗ -> (sum(λ⃗ .^ u[3]) - 3) / u[2], p.test.data.λ)
#         res .= [max_sum]
#         res
#     end
#     return (cons=f, lcons=[-Inf], ucons=[0.0])
# end

struct ValanisLandel{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = 2\\mu\\sum\\limits_{1}^{3}(\\lambda_i(\\log\\lambda_i -1))
```

# Parameters:
- `μ`

> Valanis KC, Landel RF. The strain‐energy function of a hyperelastic material in terms of the extension ratios. Journal of Applied Physics. 1967 Jun;38(7):2997-3002.
"""
ValanisLandel() = ValanisLandel{PrincipalValueForm}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::ValanisLandel{T},
    λ⃗::Vector{S},
    (; μ),
) where {T,S}
    return 2 * μ * sum(λ⃗ .* (log.(λ⃗) .- 1))
end

function parameters(::ValanisLandel)
    return (:μ,)
end

struct PengLandel{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = E\\sum\\limits_{i=1}^{3}\\bigg[\\lambda_i - 1 - \\log(\\lambda_i) - \\frac{1}{6}\\log(\\lambda_i)^2 + \\frac{1}{18}\\log(\\lambda_i)^3-\\frac{1}{216}\\log(\\lambda_i)^4\\bigg]
```

# Parameters:
- `E`

> Peng TJ, Landel RF. Stored energy function of rubberlike materials derived from simple tensile data. Journal of Applied Physics. 1972 Jul;43(7):3064-7.
"""
PengLandel() = PengLandel{PrincipalValueForm}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::PengLandel{T},
    λ⃗::Vector{S},
    (; E),
) where {T,S}
    return sum(
        @. (
            λ⃗ - 1 - log(λ⃗) - 1 / 6 * log(λ⃗)^2 + 1 / 18 * log(λ⃗)^3 - 1 / 216 * log(λ⃗)^4
        ) * E
    )
end

function parameters(::PengLandel)
    return (:E,)
end

struct Ogden{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = \\sum\\limits_{i=1}^{N}\\frac{\\mu_i}{\\alpha_i}(\\lambda_1^{\\alpha_i}+\\lambda_2^{\\alpha_i}+\\lambda_3^{\\alpha_i}-3)
```

# Parameters:
- μ⃗
- α⃗

> Ogden RW. Large deformation isotropic elasticity–on the correlation of theory and experiment for incompressible rubberlike solids. Proceedings of the Royal Society of London. A. Mathematical and Physical Sciences. 1972 Feb 1;326(1567):565-84.
"""
Ogden() = Ogden{PrincipalValueForm}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Ogden{T},
    λ⃗::Vector{S},
    (; μ⃗, α⃗),
) where {T,S}
    λ_a = λ⃗[1] .^ α⃗ + λ⃗[2] .^ α⃗ + λ⃗[3] .^ α⃗
    return sum(@. μ⃗ / α⃗ * (λ_a - 3))
end

function parameters(::Ogden)
    return (:μ⃗, :α⃗)
end

struct Attard{T} <: AbstractIncompressibleModel{T}
    Wi::Function
end

"""
$(SIGNATURES)

# Model:

```math
W = \\sum\\limits_{i=1}^N\\frac{A_i}{2i}(\\lambda_1^{2i}+\\lambda_2^{2i}+\\lambda_3^{2i}-3) + \\frac{B_i}{2i}(\\lambda_1^{-2i}+\\lambda_2^{-2i}+\\lambda_3^{-2i}-3)
```

# Parameters:
- `A⃗`
- `B⃗`

> Attard MM, Hunt GW. Hyperelastic constitutive modeling under finite strain. International Journal of Solids and Structures. 2004 Sep 1;41(18-19):5327-50.
"""
function Attard()
    f(i, (; λ⃗, p)) =
        p.A⃗[i] / 2 / i * (sum(λ⃗ .^ (2i)) - 3) + p.B⃗[i] / 2 / i * (sum(λ⃗ .^ (-2i)) - 3)
    Attard{PrincipalValueForm}(f)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::Attard{T},
    λ⃗::Vector{S},
    p,
) where {T,S}
    @assert length(p.A⃗) == length(p.B⃗) "Length of A and B are not equal"
    W = sum(Base.Fix2(ψ.Wi, (λ⃗ = λ⃗, p = p)), 1:length(p.A⃗))
    return W
end

function parameters(::Attard)
    return (:A⃗, :B⃗)
end

struct Shariff{T} <: AbstractIncompressibleModel{T}
    ϕ::Vector{Function}
    Φ::Vector{Function}
end

"""
$(SIGNATURES)

Model:

```math
W = E\\sum\\limits_{i=1}^3\\sum\\limits_{j=1}^{N}|\\alpha_j| \\Phi_j(\\lambda_i)
```

Parameters:
- `E`
- `α⃗`

> Shariff MH. Strain energy function for filled and unfilled rubberlike material. Rubber chemistry and technology. 2000 Mar;73(1):1-8.
"""
function Shariff()
    ϕ1(x) = 2 * log(x) / 3
    ϕ2(x) = exp(1 - x) + x - 2
    ϕ3(x) = exp(x - 1) - x
    ϕ4(x) = (x - 1)^3 / x^3.6
    ϕj(x, j) = (x - 1)^(j - 1)

    ϕ = [ϕ1, ϕ2, ϕ3, ϕ4, ϕj]

    c(j, r) = factorial(j) / factorial(r) / factorial(j - r)
    Φ1(x) = log(x)^2 / 3

    Φ2(x) = -exp(1.0) * expinti(-1.0) + exp(1.0) * expinti(-x) + x - 2 * log(x) - 1

    Φ3(x) = (expinti(x) - expinti(1.0)) / exp(1.0) - x + 1

    # # Φ4(x) = -1 / (0.6 * x^(0.6)) + 3 / (1.6 * x^(1.6)) - 3 / (2.6 * x^(2.6)) + 1 / (5.6 * x^(5.6)) + 107200 / 139776
    Φ4(x) = 5 / 936 * (125 + (52 - 216 * x + 351 * (x^2) - 312 * (x^3)) / (x^(18 / 5)))

    Φj(x, j) =
        (-1)^(j - 1) * log(x) +
        (-1)^(j - 1) * sum(r -> (-1)^r * c(j - 1, r) * x^r / r, range(1, j - 1)) -
        (-1)^(j - 1) * sum(r -> (-1)^r * c(j - 1, r) / r, range(1, j - 1))

    Φ = [Φ1, Φ2, Φ3, Φ4, Φj]

    return Shariff{PrincipalValueForm}(ϕ, Φ)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::Shariff{T},
    λ⃗::Vector{S},
    (; E, α⃗),
) where {T<:PrincipalValueForm,S}
    n = length(α⃗)
    W1 = sum(i -> sum(α⃗[i] * ψ.Φ[i].(λ⃗)), 1:minimum([4, n]))
    W2 = sum(i -> sum(α⃗[i] * ψ.Φ[5].(λ⃗, i)), minimum([5, n]):n)
    W = W1 + W2
    return E * W
end

function ContinuumMechanicsBase.CauchyStressTensor(
    ψ::Shariff{T},
    λ⃗::Vector{S},
    (; E, α⃗);
    kwargs...,
) where {T<:PrincipalValueForm,S}
    n = length(α⃗)
    σ1 = sum(i -> α⃗[i] .* ψ.ϕ[i].(λ⃗), 1:minimum([4, n]))
    σ2 = sum(i -> α⃗[i] .* ψ.ϕ[5].(λ⃗, i), minimum([5, n]):n)
    σ = σ1 + σ2
    return E .* σ
end

function ContinuumMechanicsBase.SecondPiolaKirchoffStressTensor(
    ψ::Shariff{T},
    λ⃗::Vector{S},
    (; E, α⃗);
    kwargs...,
) where {T<:PrincipalValueForm,S}
    n = length(α⃗)
    s1 = sum(i -> α⃗[i] .* ψ.ϕ[i].(λ⃗), 1:minimum([4, n]))
    s2 = sum(i -> α⃗[i] .* ψ.ϕ[5].(λ⃗, i), minimum([5, n]):n)
    s = s1 + s2
    return E .* s ./ λ⃗
end

function parameters(::Shariff)
    return (:E, :α⃗)
end

struct ArmanNarooei{T} <: AbstractIncompressibleModel{T}
    Wi::Function
end

"""
$(SIGNATURES)

# Model:

```math
W = \\sum\\limits_{i=1}^{N} A_i\\big[\\exp{m_i(\\lambda_1^{\\alpha_i}+\\lambda_2^{\\alpha_i}+\\lambda_3^{\\alpha_i}-3)}-1] + B_i\\big[\\exp{n_i(\\lambda_1^{-\\beta_i}+\\lambda_2^{-\\beta_i}+\\lambda_3^{-\\beta_i}-3)}-1]
```

# Parameters:
- `A⃗`
- `B⃗`
- `m⃗`
- `n⃗`
- `α⃗`
- `β⃗`

> Narooei K, Arman M. Modification of exponential based hyperelastic strain energy to consider free stress initial configuration and Constitutive modeling. Journal of Computational Applied Mechanics. 2018 Jun 1;49(1):189-96.
"""
function ArmanNarooei()
    f(i, (; λ⃗, p)) =
        p.A⃗[i] * (exp(p.m⃗[i] * (sum(λ⃗ .^ p.α⃗[i]) - 3)) - 1) +
        p.B⃗[i] * (exp(p.n⃗[i] * (sum(λ⃗ .^ (-p.β⃗[i])) - 3)) - 1)
    ArmanNarooei{PrincipalValueForm}(f)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::ArmanNarooei{T},
    λ⃗::Vector{S},
    p,
) where {T,S}
    @assert length(p.A⃗) ==
            length(p.B⃗) ==
            length(p.m⃗) ==
            length(p.n⃗) ==
            length(p.α⃗) ==
            length(p.β⃗) "Length of A, B, m, n, α. and β are not equal"
    # (; A⃗, B⃗, m⃗, n⃗, α⃗, β⃗)
    W = sum(Base.Fix2(ψ.Wi, (λ⃗ = λ⃗, p = p)), 1:length(p.A⃗))
    # @tullio W := A⃗[i] * (exp(m⃗[i] * (sum(λ⃗ .^ α⃗[i]) - 3)) - 1) + B⃗[i] * (exp(n⃗[i] * (sum(λ⃗ .^ (-β⃗[i])) - 3)) - 1)
    return W
end

function parameters(::ArmanNarooei)
    return (:A⃗, :B⃗, :m⃗, :n⃗, :α⃗, :β⃗)
end

struct ContinuumHybrid{T} <: AbstractIncompressibleModel{T} end

"""
$(SIGNATURES)

# Model:

```math
W = K_1(I_1-3)+K_2\\log\\frac{I_2}{3}+\\frac{\\mu}{\\alpha}(\\lambda_1^\\alpha+\\lambda_2^\\alpha+\\lambda^\\alpha-3)
```

# Parameters:
- `K₁`
- `K₂`
- `α`
- `μ`

> Beda T, Chevalier Y. Hybrid continuum model for large elastic deformation of rubber. Journal of applied physics. 2003 Aug 15;94(4):2701-6.
"""
ContinuumHybrid() = ContinuumHybrid{PrincipalValueForm}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::ContinuumHybrid{T},
    λ⃗::Vector{S},
    (; K₁, K₂, α, μ),
) where {T,S}
    return K₁ * (I₁(λ⃗) - 3) + K₂ * log(I₂(λ⃗) / 3) + μ / α * (sum(λ⃗ .^ α) - 3)
end

function parameters(::ContinuumHybrid)
    return (:K₁, :K₂, :α, :μ)
end

struct Bechir4Term{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = C_1^1(I_1-3)+\\sum\\limits_{n=1}^{2}\\sum\\limits_{r=1}^{2}C_n^{r}(\\lambda_1^{2n}+\\lambda_2^{2n}+\\lambda_3^{2n}-3)^r
```

# Parameters:
- `C11`
- `C12`
- `C21`
- `C22`

> Khajehsaeid H, Arghavani J, Naghdabadi R. A hyperelastic constitutive model for rubber-like materials. European Journal of Mechanics-A/Solids. 2013 Mar 1;38:144-51.
"""
Bechir4Term() = Bechir4Term{PrincipalValueForm}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Bechir4Term{T},
    λ⃗::Vector{S},
    (; C11, C12, C21, C22),
) where {T,S}
    C = [C11 C12; C21 C22]
    W1 = C[1, 1] * (I₁(λ⃗) - 3)
    W2 =
        C[1, 1] * (sum(λ⃗ .^ (2 * 1))) +
        C[1, 2] * (sum(λ⃗ .^ (2 * 1))) +
        C[2, 1] * (sum(λ⃗ .^ (2 * 2))) +
        C[2, 2] * (sum(λ⃗ .^ (2 * 2)))
    # return C[1, 1] * (I₁(λ⃗) - 3) + sum(n -> sum(r -> C[n, r] * (sum(λ⃗ .^ (2n))), 1:2), 1:2)
    return W1 + W2
end

function parameters(::Bechir4Term)
    return (:C11, :C12, :C21, :C22)
end

struct ConstrainedJunction{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = G_c (I_1-3)+ \\frac{\\nu k T}{2}\\left(\\sum\\limits_{i=1}^{3}\\kappa\\frac{\\lambda_i-1}{\\lambda_i^2+\\kappa}+\\log{\\frac{\\lambda_i^2+\\kappa}{1+\\kappa}}-\\log{\\lambda_i^2}\\right)
```

# Parameters:
- `Gc`
- `νkT`
- `κ`

> Flory PJ, Erman B. Theory of elasticity of polymer networks. 3. Macromolecules. 1982 May;15(3):800-6.
> Erman B, Flory PJ. Relationships between stress, strain, and molecular constitution of polymer networks. Comparison of theory with experiments. Macromolecules. 1982 May;15(3):806-11.
"""
ConstrainedJunction() = ConstrainedJunction{PrincipalValueForm}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::ConstrainedJunction{T},
    λ⃗::Vector{S},
    (; Gc, μkT, κ),
) where {T,S}
    return Gc * (I₁(λ⃗) - 3) +
           μkT / 2 * sum(
        i -> κ * (λ⃗[i] - 1) / (λ⃗[i]^2 + κ) + log((λ⃗[i]^2 + κ) / (1 + κ)) - log(λ⃗[i]^2),
        1:3,
    )
end

function parameters(::ConstrainedJunction)
    return (:Gc, :μkT, :κ)
end

function parameter_bounds(::ConstrainedJunction, data::AbstractHyperelasticTest)
    λ_min = minimum(minimum.(collect.(data.data.λ)))
    κ_min = -λ_min^2
    lb = (Gc = -Inf, μkT = -Inf, κ = κ_min)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct EdwardVilgis{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = \\frac{1}{2}N_C\\Bigg[\\frac{(1-\\alpha^2)I_1}{1-\\alpha^2I_1}+\\log(1-\\alpha^2I_1)\\Bigg]+\\frac{1}{2}N_S\\Bigg[\\sum_{i=1}^{3}\\Big\\{\\frac{(1+\\eta)(1-\\alpha^2)\\lambda_i^2}{( 1+\\eta\\lambda_i^2)(1-\\alpha^2I_1)}+\\log(1+\\eta\\lambda_i^2)\\Big\\}+\\log(1-\\alpha^2I_1)\\Bigg]
```

# Parameters:
- `Ns`: Number of sliplinks
- `Nc`: Number of crosslinks
- `α`: A measure of chain inextensibility
- `η`: A measure of the amount of chain slippage

# Note:
- Since α and η result from the same mechanism, they should be of approximately the same order of magnitude. Large differences between the two may indicate an issue with the optimizer or initial guess.

> Edwards SF, Vilgis T. The effect of entanglements in rubber elasticity. Polymer. 1986 Apr 1;27(4):483-92.
"""
EdwardVilgis() = EdwardVilgis{PrincipalValueForm}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::EdwardVilgis{T},
    λ⃗::Vector{S},
    (; Ns, Nc, α, η),
) where {T,S}
    W_Nc = 0.5 * Nc * ((1 - α^2) * I₁(λ⃗) / (1 - α^2 * I₁(λ⃗)) + log(1 - α^2 * I₁(λ⃗)))
    W_Ns =
        0.5 *
        Ns *
        (
            (1 + η) * (1 - α^2) * λ⃗[1] / (1 + η * λ⃗[1]^2) / (1 - α^2 * I₁(λ⃗)) +
            log(1 + η * λ⃗[1]^2) +
            (1 + η) * (1 - α^2) * λ⃗[2] / (1 + η * λ⃗[2]^2) / (1 - α^2 * I₁(λ⃗)) +
            log(1 + η * λ⃗[2]^2) +
            (1 + η) * (1 - α^2) * λ⃗[3] / (1 + η * λ⃗[3]^2) / (1 - α^2 * I₁(λ⃗)) +
            log(1 + η * λ⃗[3]^2) +
            log(1 - α^2 * I₁(λ⃗))
        )
    W = W_Nc + W_Ns
    return W
end

function parameters(::EdwardVilgis)
    return (:Ns, :Nc, :α, :η)
end

function parameter_bounds(::EdwardVilgis, data::AbstractHyperelasticTest)
    # I₁_max = maximum()
    λ_max = maximum(maximum.(data.data.λ))
    η_min = -1 / λ_max^2
    α_max = minimum(@. sqrt(1 / I₁(data.data.λ)))
    lb = (Ns = -Inf, Nc = -Inf, α = 0.0, η = 0.0)
    ub = (Ns = Inf, Nc = Inf, α = α_max, η = Inf)
    return (lb = lb, ub = ub)
end

struct MCC{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = \\frac{1}{2}\\zeta k T \\sum\\limits_{i=1}^{3}(\\lambda_i^2-1)+\\frac{1}{2}\\mu k T\\sum\\limits_{i=1}^{3}[B_i+D_i-\\log{(1+B_i)}-\\log{(1+D_i)}]
```

where:

```math
B_i = \\frac{\\kappa^2(\\lambda_i^2-1)}{(\\lambda_i^2+\\kappa)^2}
```

and

```math
D_i = \\frac{\\lambda_i^2 B_i}{\\kappa}
```

# Parameters:
- `ζkT`
- `μkT`
- `κ`

> Erman B, Monnerie L. Theory of elasticity of amorphous networks: effect of constraints along chains. Macromolecules. 1989 Aug;22(8):3342-8.
"""
MCC() = MCC{PrincipalValueForm}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::MCC{T},
    λ⃗::Vector{S},
    (; ζkT, μkT, κ),
) where {T,S}
    B = @. κ^2 * (λ⃗^2 - 1) * (λ⃗^2 + κ)^(-2)
    D = @. λ⃗^2 * B / κ
    W1 = @. λ⃗^2 - 1
    W2 = @. B - log(1 + B)
    W3 = @. D - log(1 + D)
    return sum(1 / 2 * ζkT * W1 + 1 / 2 * μkT * (W2 + W3))
end

function parameters(::MCC)
    return (:ζkT, :μkT, :κ)
end

function parameter_bounds(::MCC, data::AbstractHyperelasticTest)
    lb = (ζkT = -Inf, μkT = -Inf, κ = 0)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct Tube{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = \\sum\\limits_{i=1}^{3}\\frac{G_c}{2}(\\lambda_i^2-1)+\\frac{2Ge}{\\beta^2}(\\lambda_i^{-\\beta}-1)
```

# Parameters:
- `Gc`
- `Ge`
- `β`

> Heinrich G, Kaliske M. Theoretical and numerical formulation of a molecular based constitutive tube-model of rubber elasticity. Computational and Theoretical Polymer Science. 1997 Jan 1;7(3-4):227-41.
"""
Tube() = Tube{PrincipalValueForm}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::Tube{T},
    λ⃗::Vector{S},
    (; Gc, Ge, β),
) where {T,S}
    return sum(@. Gc / 2 * (λ⃗^2 - 1) + 2Ge / β^2 * (λ⃗^(-β) - 1))
end

function parameters(::Tube)
    return (:Gc, :Ge, :β)
end

struct NonaffineTube{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = G_c \\sum\\limits_{i=1}^{3}\\frac{\\lambda_i^2}{2}+G_e\\sum\\limits_{i=1}^{3}\\lambda_i+\\frac{1}{\\lambda_i}
```

# Parameters:
- `Gc`
- `Ge`

> Rubinstein M, Panyukov S. Nonaffine deformation and elasticity of polymer networks. Macromolecules. 1997 Dec 15;30(25):8036-44.
"""
NonaffineTube() = NonaffineTube{PrincipalValueForm}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::NonaffineTube{T},
    λ⃗::Vector{S},
    (; Gc, Ge),
) where {T,S}
    return Gc * sum(λ⃗ .^ 2 ./ 2) + Ge * sum(λ⃗ .+ 1 ./ λ⃗)
end

function parameters(::NonaffineTube)
    return (:Gc, :Ge)
end

struct ThreeChainModel{T} <: AbstractIncompressibleModel{T}
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation
end

"""
$(SIGNATURES)

# Model:

```math
W = \\frac{\\mu\\sqrt{N}}{3}\\sum\\limits_{i=1}^{3}\\bigg(\\lambda_i\\beta_i+\\sqrt{N}\\log\\bigg(\\frac{\\beta_i}{\\sinh \\beta_i}\\bigg)\\bigg)
```

# Arguments:
- `ℒinv=TreloarApproximation()`: Sets the inverse Langevin approxamation used

# Parameters:
- `μ`: Small strain shear modulus
- `N`: Square of the locking stretch of the network.

> James HM, Guth E. Theory of the elastic properties of rubber. The Journal of Chemical Physics. 1943 Oct;11(10):455-81.
"""
ThreeChainModel(;
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation = TreloarApproximation(),
) = ThreeChainModel{PrincipalValueForm}(ℒinv)

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::ThreeChainModel{T},
    λ⃗::Vector{S},
    (; μ, N),
) where {T,S}
    β = @. inverse_langevin_approximation(λ⃗ / sqrt(N), ψ.ℒinv)
    return μ * sqrt(N) / 3 * sum(@. λ⃗ * β + sqrt(N) * log(β / sinh(β)))
end

function parameters(::ThreeChainModel)
    return (:μ, :N)
end

function parameter_bounds(::ThreeChainModel, data::AbstractHyperelasticTest)
    λ_max = maximum(maximum.(collect.(data.data.λ)))
    N_min = λ_max^2
    lb = (μ = -Inf, N = N_min)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct ModifiedFloryErman{T} <: AbstractIncompressibleModel{T}
    # ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation
    Chain8::ArrudaBoyce
end
"""
$(SIGNATURES)

# Model:

```math
W = W_{\\text{Arruda-Boyce}}+\\sum\\limits_{i=1}^{3}\\frac{\\mu}{2}[B_i+D_i]
```

# Arguments:
- `ℒinv=TreloarApproximation()`: Sets the inverse Langevin approxamation used

# Parameters:
- `μ`
- `N`
- `κ`

> Edwards SF. The statistical mechanics of polymerized material. Proceedings of the Physical Society (1958-1967). 1967 Sep 1;92(1):9.
"""
function ModifiedFloryErman(;
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation = TreloarApproximation(),
)
    ModifiedFloryErman{PrincipalValueForm}(ArrudaBoyce(PrincipalValueForm(), ℒinv = ℒinv))
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    W::ModifiedFloryErman{T},
    λ⃗::Vector{S},
    p,
) where {T,S}
    WAB = StrainEnergyDensity(W.Chain8, λ⃗, p)
    B = @. p.κ^2 * (λ⃗^2 - 1) / (λ⃗^2 + p.κ)^2
    D = @. λ⃗^2 * B / p.κ
    W2 = sum(@. B + D - log(B + 1) - log(D + 1))
    return WAB + W2
end

function parameters(::ModifiedFloryErman)
    return (:μ, :N, :κ)
end

function parameter_bounds(::ModifiedFloryErman, data::AbstractHyperelasticTest)
    I₁_max = maximum(I₁.(data.data.λ))
    # N_max = 11 / 35 * I₁_max # old
    N_max = I₁_max / 3
    lb = (μ = -Inf, N = N_max, κ = 0)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct ExtendedTubeModel{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = \\frac{G_c}{2}\\bigg[\\frac{(1-\\delta^2)(I_1-3)}{1-\\delta^2(I_1-3)}+\\log{(1-\\delta^2(I_1-3))}\\bigg]+\\frac{2G_e}{\\beta^2}\\sum\\limits_{i=1}^{3}(\\lambda_i^{-\\beta}-1)
```

# Parameters:
- `Gc`
- `Ge`
- `δ`
- `β`

> Kaliske M, Heinrich G. An extended tube-model for rubber elasticity: statistical-mechanical theory and finite element implementation. Rubber Chemistry and Technology. 1999 Sep;72(4):602-32.
"""
ExtendedTubeModel() = ExtendedTubeModel{PrincipalValueForm}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::ExtendedTubeModel{T},
    λ⃗::Vector{S},
    (; Gc, Ge, δ, β),
) where {T,S}
    return Gc / 2 * (
        (1 - δ^2) * (I₁(λ⃗) - 3) / (1 - δ^2 * (I₁(λ⃗) - 3)) + log(1 - δ^2 * (I₁(λ⃗) - 3))
    ) + 2 * Ge / β^2 * sum(λ⃗ .^ (-β) .- 1)
end

function parameters(::ExtendedTubeModel)
    return (:Gc, :Ge, :δ, :β)
end

function parameter_bounds(::ExtendedTubeModel, data::AbstractHyperelasticTest)
    I₁_max = maximum(I₁.(data.data.λ))

    δ_max = sqrt(1 / (I₁_max - 3))
    lb = (Gc = -Inf, Ge = -Inf, δ = -δ_max, β = 0)
    ub = (Gc = Inf, Ge = Inf, δ = δ_max, β = Inf)
    return (lb = lb, ub = ub)
end

struct NonaffineMicroSphere{T} <: AbstractIncompressibleModel{T}
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation
    r⃗²::Vector{Vector{Float64}}
    w::Vector{Float64}
end

"""
$(SIGNATURES)

# Model: See Paper

# Arguments:
- `ℒinv=CohenRounded3_2()`: Sets the inverse Langevin approximation used.
- `n=21`: Order of quadrature for spherical integration

# Parameters:
- μ: Small strain shear modulus
- N: Number of chain segments
- p: Non-affine stretch parameter
- U: Tube geometry parameter
- q: Non-affine tube parameter

> Miehe C, Göktepe S, Lulei F. A micro-macro approach to rubber-like materials—part I: the non-affine micro-sphere model of rubber elasticity. Journal of the Mechanics and Physics of Solids. 2004 Nov 1;52(11):2617-60.
"""
function NonaffineMicroSphere(;
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation = CohenRounded3_2(),
    n = 21,
)
    a = √(2) / 2
    b = 0.836095596749
    c = 0.387907304067
    r⃗² = [
        [0, 0, 1] .^ 2,
        [0, 1, 0] .^ 2,
        [1, 0, 0] .^ 2,
        [0, a, a] .^ 2,
        [0, -a, a] .^ 2,
        [a, 0, a] .^ 2,
        [-a, 0, a] .^ 2,
        [a, a, 0] .^ 2,
        [-a, a, 0] .^ 2,
        [b, c, c] .^ 2,
        [-b, c, c] .^ 2,
        [b, -c, c] .^ 2,
        [-b, -c, c] .^ 2,
        [c, b, c] .^ 2,
        [-c, b, c] .^ 2,
        [c, -b, c] .^ 2,
        [-c, -b, c] .^ 2,
        [c, c, b] .^ 2,
        [-c, c, b] .^ 2,
        [c, -c, b] .^ 2,
        [-c, -c, b] .^ 2,
    ]
    w1 = 0.02652142440932
    w2 = 0.0199301476312
    w3 = 0.0250712367487
    w = 2 .* [fill(w1, 3); fill(w2, 6); fill(w3, 12)] # Multiply by two since integration is over the half-sphere
    NonaffineMicroSphere{PrincipalValueForm}(ℒinv, r⃗², w)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ψ::NonaffineMicroSphere{T},
    λ⃗::Vector{S},
    (; μ, N, p, U, q),
) where {T,S}
    #
    λ⃗² = λ⃗ .^ 2
    inv_λ⃗² = inv.(λ⃗²)

    λ⃗²_r⃗² = broadcast(
        Base.Fix2((x, y) -> sqrt(x[1] * y[1] + x[2] * y[2] + x[3] * y[3]), λ⃗²),
        ψ.r⃗²,
    )

    inv_λ⃗²_r⃗² = broadcast(
        Base.Fix2((x, y) -> sqrt(x[1] * y[1] + x[2] * y[2] + x[3] * y[3]), inv_λ⃗²),
        ψ.r⃗²,
    )

    λ = sum((λ⃗²_r⃗² .^ p) .* ψ.w)
    λr = λ^(1 / p) / √N
    β = inverse_langevin_approximation(λr, ψ.ℒinv)
    ψf = μ * N * (λr * β + log(β / sinh(β)))

    ν = sum((sqrt.(inv_λ⃗²_r⃗²) .^ q) .* ψ.w)
    ψc = U * μ * N * (ν)

    return ψf + ψc
end

function parameters(::NonaffineMicroSphere)
    return (:μ, :N, :p, :U, :q)
end

function parameter_bounds(::NonaffineMicroSphere, data::AbstractHyperelasticTest)
    lb = (μ = -Inf, N = 0, p = 0, U = 0, q = 0)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct Bootstrapped8Chain{T} <: AbstractIncompressibleModel{T}
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation
    W8::Function
end

"""
$(SIGNATURES)

# Model:

```math
W = W_8(\\frac{\\sum\\lambda}{\\sqrt{3N}}-\\frac{\\lambda_{chain}}{\\sqrt{N}})+W_{8}(\\frac{\\lambda_{chain}}{\\sqrt{N}})
```

where:

```math
W_8(x) = \\mu N (x \\mathcal{L}^{-1}(x) + \\log\\frac{\\mathcal{L}^{-1}(x)}{\\sinh\\mathcal{L}^{-1}(x)})
```

and

```math
\\lambda_{chain} = \\sqrt{\\frac{I_1}{3}}
```

# Arguments:
- `ℒinv=TreloarApproximation()`: Sets the inverse Langevin approximation used.

# Parameters:
- `μ`
- `N`

> Miroshnychenko D, Green WA, Turner DM. Composite and filament models for the mechanical behaviour of elastomeric materials. Journal of the Mechanics and Physics of Solids. 2005 Apr 1;53(4):748-70.
> Miroshnychenko D, Green WA. Heuristic search for a predictive strain-energy function in nonlinear elasticity. International Journal of Solids and Structures. 2009 Jan 15;46(2):271-86.

"""
function Bootstrapped8Chain(;
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation = TreloarApproximation(),
)
    function W8(x, (; μ, N))
        β = inverse_langevin_approximation(x, ℒinv)
        μ * N * (x * β + log(β / sinh(β)))
    end
    Bootstrapped8Chain{PrincipalValueForm}(ℒinv, W8)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    W::Bootstrapped8Chain,
    λ⃗::Vector{T},
    p,
) where {T}
    λchain = √(I₁(λ⃗) / 3)
    W.W8(sum(λ⃗) / √(3 * p.N) - λchain / √(p.N), p) + W.W8(λchain / √(p.N), p)
end

function parameters(::Bootstrapped8Chain)
    return (:μ, :N)
end

function parameter_bounds(::Bootstrapped8Chain, data::AbstractHyperelasticTest)
    I₁_max = maximum(I₁.(data.data.λ))
    N_min = I₁_max / 3
    lb = (μ = -Inf, N = N_min)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct DavidsonGoulbourne{T} <: AbstractIncompressibleModel{T} end

"""
$(SIGNATURES)

# Model:

```math
W = \\frac{G_c I_1}{6}-G_c\\lambda_{max}\\log\\left(3\\lambda_{max}^2-I_1\\right)+G_e\\sum\\limits_{i=1}^{3}\\left(\\lambda_i+\\frac{1}{\\lambda_i}\\right)
```

# Parameters:
- Gc
- Ge
- λmax

> Davidson JD, Goulbourne NC. A nonaffine network model for elastomers undergoing finite deformations. Journal of the Mechanics and Physics of Solids. 2013 Aug 1;61(8):1784-97.
"""
DavidsonGoulbourne() = DavidsonGoulbourne{PrincipalValueForm}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::DavidsonGoulbourne{T},
    λ⃗::Vector{S},
    (; Gc, Ge, λmax),
) where {T,S}
    return 1 / 6 * Gc * I₁(λ⃗) - Gc * λmax^2 * log(3 * λmax^2 - I₁(λ⃗)) +
           Ge * (λ⃗[1] + 1 / λ⃗[1] + λ⃗[2] + 1 / λ⃗[2] + λ⃗[3] + 1 / λ⃗[3])
end

function parameters(::DavidsonGoulbourne)
    return (:Gc, :Ge, :λmax)
end

function parameter_bounds(::DavidsonGoulbourne, data::AbstractHyperelasticTest)
    I₁_max = maximum(I₁.(data.data.λ))
    λmax_min = sqrt(I₁_max / 3)
    lb = (Gc = 0, Ge = 0, λmax = λmax_min)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct KhiemItskov{T} <: AbstractIncompressibleModel{T} end
"""
$(SIGNATURES)

# Model:

```math
W = \\mu_c \\kappa n \\log\\bigg(\\frac{\\sin(\\frac{\\pi}{\\sqrt{n}})(\\frac{I_1}{3})^{\\frac{q}{2}}}{\\sin(\\frac{\\pi}{\\sqrt{n}}(\\frac{I_1}{3})^{\\frac{q}{2}}}\\bigg)+\\mu_t\\big[\\frac{I_2}{3}^{1/2} - 1 \\big]
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()

# Parameters:
- μcκ
- n
- q
- μt

> Khiêm VN, Itskov M. Analytical network-averaging of the tube model:: Rubber elasticity. Journal of the Mechanics and Physics of Solids. 2016 Oct 1;95:254-69.
"""
KhiemItskov(
    type::T = PrincipalValueForm(),
) where {T<:Union{PrincipalValueForm,InvariantForm}} = KhiemItskov{T}()

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::KhiemItskov{T},
    λ⃗::Vector{S},
    (; μcκ, n, q, μt),
) where {T<:PrincipalValueForm,S}
    I1 = I₁(λ⃗)
    num = (sin(π / sqrt(n)) * (I1 / 3)^(q / 2))
    denom = (sin(π / sqrt(n) * (I1 / 3)^(q / 2)))
    # @assert num ≥ denom "Parameters are not feasible"
    return μcκ * n * log(num / denom) + μt * ((I₂(λ⃗) / 3)^(1 / 2) - 1)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::KhiemItskov{T},
    I⃗::Vector{S},
    (; μcκ, n, q, μt),
) where {T<:InvariantForm,S}
    num = (sin(π / sqrt(n)) * (I⃗[1] / 3)^(q / 2))
    denom = (sin(π / sqrt(n) * (I⃗[1] / 3)^(q / 2)))
    # @assert num ≥ denom "Parameters are not feasible: $((μcκ, n, q, μt))"
    return μcκ * n * log(num / denom) + μt * ((I⃗[2] / 3)^(1 / 2) - 1)
end

function parameters(::KhiemItskov)
    return (:μcκ, :n, :q, :μt)
end

function parameter_bounds(::KhiemItskov, data::AbstractHyperelasticTest)
    lb = (n = 0, μcκ = -Inf, μt = Inf, q = 0)
    ub = (n = Inf, μcκ = -Inf, μt = Inf, q = Inf)
    return (lb = lb, ub = ub)
end
# function constraints(::KhiemItskov, data::AbstractHyperelasticTest)
#     I₁_max = maximum(I₁.(data.data.λ))
#     f(u, p) = [(sin(π / sqrt(u.n)) * (I₁_max / 3)^(u.q / 2)) / (sin(π / sqrt(u.n) * (I₁_max / 3)^(u.q / 2)))]
#     return f
# end


struct GeneralConstitutiveModel_Network{T} <: AbstractIncompressibleModel{T}
    GeneralConstitutiveModel_Network(
        ::T = PrincipalValueForm(),
    ) where {T<:Union{PrincipalValueForm,InvariantForm}} = new{T}()
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::GeneralConstitutiveModel_Network{T},
    λ⃗::Vector{S},
    (; Gc, N),
) where {T<:PrincipalValueForm,S}
    I1 = I₁(λ⃗)
    return Gc * N * log((3 * N + 0.5 * I1) / (3 * N - I1))
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::GeneralConstitutiveModel_Network{T},
    I⃗::Vector{S},
    (; Gc, N),
) where {T<:InvariantForm,S}
    return Gc * N * log((3 * N + 0.5 * I⃗[1]) / (3 * N - I⃗[1]))
end


function parameters(::GeneralConstitutiveModel_Network)
    return (:Gc, :N)
end

function parameter_bounds(
    ::GeneralConstitutiveModel_Network,
    data::AbstractHyperelasticTest,
)
    I₁_max = maximum(I₁.(data.data.λ))
    N_min = I₁_max / 3
    lb = (Gc = -Inf, N = N_min)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct GeneralConstitutiveModel_Tube{T} <: AbstractIncompressibleModel{T}
    GeneralConstitutiveModel_Tube(
        ::T = PrincipalValueForm(),
    ) where {T<:PrincipalValueForm} = new{T}()
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    ::GeneralConstitutiveModel_Tube{T},
    λ⃗::Vector{S},
    (; Ge),
) where {T,S}
    return sum(Ge ./ λ⃗)
end

function parameters(::GeneralConstitutiveModel_Tube)
    return (:Ge,)
end

struct GeneralConstitutiveModel{T} <: AbstractIncompressibleModel{T}
    Tube::GeneralConstitutiveModel_Tube
    Network::GeneralConstitutiveModel_Network
end

"""
$(SIGNATURES)

# Model:

```math
W = G_c N \\log\\bigg(\\frac{3N+\\frac{1}{2}I_1}{3N-I_1}\\bigg)+G_e\\sum\\limits_{i=1}^{3}\\frac{1}{\\lambda_I}
```

# Parameters:
- `Gc`
- `Ge`
- `N`

> Xiang Y, Zhong D, Wang P, Mao G, Yu H, Qu S. A general constitutive model of soft elastomers. Journal of the Mechanics and Physics of Solids. 2018 Aug 1;117:110-22.
"""
GeneralConstitutiveModel() = GeneralConstitutiveModel{PrincipalValueForm}(
    GeneralConstitutiveModel_Tube(PrincipalValueForm()),
    GeneralConstitutiveModel_Network(PrincipalValueForm()),
)

function ContinuumMechanicsBase.StrainEnergyDensity(
    W::GeneralConstitutiveModel{T},
    λ⃗::Vector{S},
    ps,
) where {T,S}
    return StrainEnergyDensity(W.Network, λ⃗, ps) + StrainEnergyDensity(W.Tube, λ⃗, ps)
end

function parameters(W::GeneralConstitutiveModel)
    return (parameters(W.Network)..., parameters(W.Tube)...)
end

function parameter_bounds(W::GeneralConstitutiveModel, data::AbstractHyperelasticTest)
    I₁_max = maximum(I₁.(data.data.λ))
    N_min = I₁_max / 3
    lb = (Gc = -Inf, Ge = -Inf, N = N_min)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct FullNetwork{T} <: AbstractIncompressibleModel{T}
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation
    Chain3::ThreeChainModel
    Chain8::ArrudaBoyce
end

"""
$(SIGNATURES)

# Model:

```math
W = (1-\\rho)W_{3Chain}+\\rho W_{8chain}
```

# Arguments:
- `ℒinv=TreloarApproximation()`: Sets the inverse Langevin approxamation used

# Parameters:
- `μ`
- `N`
- `ρ`

> Treloar LR, Riding G. A non-Gaussian theory for rubber in biaxial strain. I. Mechanical properties. Proceedings of the Royal Society of London. A. Mathematical and Physical Sciences. 1979 Dec 31;369(1737):261-80.
> Wu PD, van der Giessen E. On improved 3-D non-Gaussian network models for rubber elasticity. Mechanics research communications. 1992 Sep 1;19(5):427-33.
> Wu PD, Van Der Giessen E. On improved network models for rubber elasticity and their applications to orientation hardening in glassy polymers. Journal of the Mechanics and Physics of Solids. 1993 Mar 1;41(3):427-56.
"""
function FullNetwork(;
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation = TreloarApproximation(),
)
    FullNetwork{PrincipalValueForm}(
        ℒinv,
        ThreeChainModel{PrincipalValueForm}(ℒinv),
        ArrudaBoyce{PrincipalValueForm}(ℒinv),
    )
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    W::FullNetwork,
    λ⃗::Vector{T},
    p,
) where {T}
    W3 = StrainEnergyDensity(W.Chain3, λ⃗, p)
    W8 = StrainEnergyDensity(W.Chain8, λ⃗, p)
    return (1 - p.ρ) * W3 + p.ρ * W8
end

function parameters(::FullNetwork)
    return (:μ, :N, :ρ)
end

function parameter_bounds(::FullNetwork, data::AbstractHyperelasticTest)
    I₁_max = maximum(I₁.(data.data.λ))
    λ_max = maximum(maximum.(data.data.λ))
    N₁ = λ_max^2
    N₂ = I₁_max / 3
    N_min = (N₁ > N₂) ? N₁ : N₂
    lb = (μ = -Inf, N = N_min, ρ = 0.0)
    ub = (μ = Inf, N = Inf, ρ = 1.0)
    return (lb = lb, ub = ub)
end

struct ZunigaBeatty{T} <: AbstractIncompressibleModel{T}
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation
    Chain3::ThreeChainModel
    Chain8::ArrudaBoyce
end

"""
$(SIGNATURES)

# Model:

```math
W = \\sqrt{\\frac{N_3+N_8}{2N_3}}W_{3Chain}+\\sqrt{\\frac{I_1}{3N_8}}W_{8Chain}
```

# Arguments:
- `ℒinv=TreloarApproximation()`: Sets the inverse Langevin approxamation used

# Parameters:
- `μ`
- `N₃`
- `N₈`

> Elı́as-Zúñiga A, Beatty MF. Constitutive equations for amended non-Gaussian network models of rubber elasticity. International journal of engineering science. 2002 Dec 1;40(20):2265-94.
"""
function ZunigaBeatty(;
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation = TreloarApproximation(),
)
    ZunigaBeatty{PrincipalValueForm}(
        ℒinv,
        ThreeChainModel{PrincipalValueForm}(ℒinv),
        ArrudaBoyce{PrincipalValueForm}(ℒinv),
    )
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    W::ZunigaBeatty{T},
    λ⃗::Vector{S},
    (; μ, N₃, N₈),
) where {T,S}
    ΛL = √((N₃ + N₈) / 2)
    ρ₃ = ΛL / √(N₃)
    W3 = StrainEnergyDensity(W.Chain3, λ⃗, (μ = μ, N = N₃))
    W8 = StrainEnergyDensity(W.Chain8, λ⃗, (μ = μ, N = N₈))
    Λch = 1 / √(3) * √(I₁(λ⃗))
    ρ₈ = Λch / √(N₈)
    return ρ₃ * W3 + ρ₈ * W8
end

function parameters(::ZunigaBeatty)
    return (:μ, :N₃, :N₈)
end

function parameter_bounds(::ZunigaBeatty, data::AbstractHyperelasticTest)
    λ_max = maximum(maximum.(data.data.λ))
    I₁_max = maximum(I₁.(data.data.λ))
    N₃_min = λ_max^2
    N₈_min = I₁_max / 3
    lb = (μ = -Inf, N₃ = N₃_min, N₈ = N₈_min)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct Lim{T} <: AbstractIncompressibleModel{T}
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation
    F::Function
    NH::NeoHookean
    AB::ArrudaBoyce
end

"""
$(SIGNATURES)

# Model:

```math
W = \\left(1-f\\left(\\frac{I_1-3}{\\hat{I_1}-3}\\right)\\right)W_{NeoHookean}(μ₁)+fW_{ArrudaBoyce}(μ₂, N)
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()
- `ℒinv=TreloarApproximation()`: Sets the inverse Langevin approxamation used

# Parameters:
- `μ₁`
- `μ₂`
- `N`
- `Î₁`

> Lim GT. Scratch behavior of polymers. Texas A&M University; 2005.
"""
function Lim(
    type::T = PrincipalValueForm();
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation = TreloarApproximation(),
) where {T<:Union{InvariantForm,PrincipalValueForm}}
    f(x) = x^3 * (10 - 15x + 6x^2)
    Lim{T}(ℒinv, f, NeoHookean(type), ArrudaBoyce(type; ℒinv))
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    W::Lim{T},
    λ⃗::Vector{S},
    (; μ₁, μ₂, N, Î₁),
) where {T<:PrincipalValueForm,S}
    Wg = StrainEnergyDensity(W.NH, λ⃗, (μ = μ₁,))
    W8 = StrainEnergyDensity(W.AB, λ⃗, (μ = μ₂, N = N))
    ζ = (I₁(λ⃗) - 3) / (Î₁ - 3)
    return (1 - W.F(ζ)) * Wg + W.F(ζ) * W8
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    W::Lim{T},
    I⃗::Vector{S},
    (; μ₁, μ₂, N, Î₁),
) where {T<:InvariantForm,S}
    Wg = StrainEnergyDensity(W.NH, I⃗, (μ = μ₁,))
    W8 = StrainEnergyDensity(W.AB, I⃗, (μ = μ₂, N = N))
    ζ = (I⃗[1] - 3) / (Î₁ - 3)
    return (1 - W.F(ζ)) * Wg + W.F(ζ) * W8
end

function parameters(::Lim)
    return (:μ₁, :μ₂, :N, :Î₁)
end

function parameter_bounds(::Lim, data::AbstractHyperelasticTest)
    I₁_max = maximum(I₁.(data.data.λ))
    N_min = I₁_max / 3
    lb = (μ₁ = -Inf, μ₂ = -Inf, N = N_min, Î₁ = 3)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct BechirChevalier{T} <: AbstractIncompressibleModel{T}
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation
    Chain3::ThreeChainModel
    Chain8::ArrudaBoyce
end

"""
$(SIGNATURES)

# Model:

```math
W = W_{3Chain}(\\mu_f, N_3)+W_{8Chain}(\\frac{\\mu_c}{3}, N_8)
```

where:

```math
\\mu_f = \\rho\\sqrt{\\frac{I_1}{3N_8}}
```

```math
\\mu_c = \\bigg(1-\\frac{\\eta\\alpha}{\\sqrt{N_3}}\\bigg)\\mu_0
```

```math
\\alpha = \\max{\\lambda_1, \\lambda_2, \\lambda_3}
```

# Arguments:
- `ℒinv=TreloarApproximation()`: Sets the inverse Langevin approxamation used


# Parameters:
- `μ₀`
- `η`
- `ρ`
- `N₃`
- `N₈`

> Bechir H, Chevalier L, Idjeri M. A three-dimensional network model for rubber elasticity: The effect of local entanglements constraints. International journal of engineering science. 2010 Mar 1;48(3):265-74.
"""
function BechirChevalier(;
    ℒinv::InverseLangevinApproximations.AbstractInverseLangevinApproximation = TreloarApproximation(),
)
    BechirChevalier{PrincipalValueForm}(
        ℒinv,
        ThreeChainModel{PrincipalValueForm}(ℒinv),
        ArrudaBoyce{PrincipalValueForm}(ℒinv),
    )
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    W::BechirChevalier{T},
    λ⃗::Vector{S},
    (; μ₀, η, ρ, N₃, N₈),
) where {T,S}
    μf = ρ * √(I₁(λ⃗) / 3 / N₈)
    W3 = StrainEnergyDensity(W.Chain3, λ⃗, (μ = μf, N = N₃))
    α = maximum(λ⃗)
    μc = (1 - η * α / √(N₃)) * μ₀
    W8 = StrainEnergyDensity(W.Chain8, λ⃗, (μ = μc / 3, N = N₈))
    return W3 + W8
end

function parameters(::BechirChevalier)
    return (:μ₀, :η, :ρ, :N₃, :N₈)
end

function parameter_bounds(::BechirChevalier, data::AbstractHyperelasticTest)
    lb = (μ₀ = -Inf, η = -Inf, ρ = -Inf, N₃ = 0, N₈ = 0)
    ub = nothing
    return (lb = lb, ub = ub)
end

struct AnsarriBenam{T} <: AbstractIncompressibleModel{T}
    n::Int
end

"""
$(SIGNATURES)

# Model:

```math
W = \\frac{3(n-1)}{2n}\\mu N \\left[\\frac{1}{3N(n-1)}(I_1 - 3) - \\log{\\frac{I_1 - 3N}{3 -3N}} \\right]
```

# Arguments:
- `type=PrincipalValueForm()`: Sets the form of the strain energy density function. Either `PrincipalValueForm()` or `InvariantForm()
- `ℒinv=TreloarApproximation()`: Sets the inverse Langevin approxamation used (default = ``)
- `n::Int=3`: Sets the order of the model

# Parameters:
- `μ`
- `n`
- `N`

> Anssari-Benam A. On a new class of non-Gaussian molecular-based constitutive models with limiting chain extensibility for incompressible rubber-like materials. Mathematics and Mechanics of Solids. 2021 Nov;26(11):1660-74.
"""
function AnsarriBenam(
    type::T = PrincipalValueForm();
    n::Int = 3,
) where {T<:Union{PrincipalValueForm,InvariantForm}}
    @assert n > 1
    AnsarriBenam{T}(n)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    W::AnsarriBenam{T},
    λ⃗::Vector{S},
    (; μ, N, C₂, γ),
) where {T<:PrincipalValueForm,S}
    return (3 * (W.n - 1)) / (2 * W.n) *
           μ *
           N *
           ((I₁(λ⃗) - 3) / (3 * N * (W.n - 1)) - log((I₁(λ⃗) - 3 * N) / (3 - 3 * N))) +
           C₂ * log((I₂(λ⃗) / 3)^γ)
end

function ContinuumMechanicsBase.StrainEnergyDensity(
    W::AnsarriBenam{T},
    I⃗::Vector{S},
    (; μ, N, C₂, γ),
) where {T<:InvariantForm,S}
    return (3 * (W.n - 1)) / (2 * W.n) *
           μ *
           N *
           ((I⃗[1] - 3) / (3N * (W.n - 1)) - log((I⃗[1] - 3N) / (3 - 3N))) +
           C₂ * log(I⃗[2] / 3)^γ
end

function parameters(::AnsarriBenam)
    return (:μ, :N, :C₂, :γ)
end

function parameter_bounds(::AnsarriBenam, test::AbstractHyperelasticTest)
    N_min = maximum(I₁, test.data.λ)
    lb = (μ = -Inf, N = N_min, C₂ = -Inf, γ = -Inf)
    ub = nothing
    return (lb = lb, ub = ub)
end
