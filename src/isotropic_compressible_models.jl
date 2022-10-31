export GeneralCompressible, LogarithmicCompressible
"""
Generic Compressible Model

Parameters:
- ψ
    - Incompressible model parameters (see model selected for parameter names)
- κ

Model: ``\\psi_{compressible} = \\psi_{incompressible}(\\vec{\\lambda}_{incompressible})+\\kappa(J-1)^2``

Example Implementation: `StrainEnergyDensityFunction(GeneralCompressible(NeoHookean()), λ⃗, (ψ = (μ=μ, ), κ = κ))`
"""
struct GeneralCompressible <: AbstractHyperelasticModel
    incompressible::AbstractHyperelasticModel
end

function NonlinearContinua.StrainEnergyDensity(ψ::GeneralCompressible, λ⃗, p)
    StrainEnergyDensity(ψ.incompressible, λ⃗, p.ψ) + p.κ / 2 * (prod(λ⃗) - 1)^2
end

function NonlinearContinua.CauchyStressTensor(ψ::GeneralCompressible, λ⃗::Vector, p; adb=AD.ForwardDiffBackend())
    σ_dev = p.κ * (prod(λ⃗) - 1)
    σ = CauchyStressTensor(ψ.incompressible, λ⃗ ./ cbrt(prod(λ⃗)), p.ψ, adb=adb)
    return σ .+ σ_dev
end

function NonlinearContinua.CauchyStressTensor(ψ::GeneralCompressible, F::Matrix, p; adb=AD.ForwardDiffBackend())
    σ_dev = p.κ * (J(F) - 1)
    σ = CauchyStressTensor(ψ.incompressible, F ./ cbrt(J(F)), p.ψ, adb=adb)
    return σ .+ σ_dev
end

function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::GeneralCompressible, λ⃗::Vector, p; adb=AD.ForwardDiffBackend())
    s_dev = p.κ * (J(F) - 1) ./ λ⃗
    s = SecondPiolaKirchoffStressTensor(ψ.incompressible, λ⃗ ./ cbrt(J(λ⃗)), p.ψ, adb=adb)
    return s .+ s_dev
end

function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::GeneralCompressible, F::Matrix, p; adb=AD.ForwardDiffBackend())
    s_dev = p.κ * (J(F) - 1) .* inv(F)
    s = SecondPiolaKirchoffStressTensor(ψ.incompressible, F ./ cbrt(J(F)), p.ψ, adb=adb)
    return s .+ s_dev
end

export GeneralCompressible
"""
Logarithmic Compressible Model

Parameters: ψ, κ

Model: ``\\psi_{compressible} = \\psi_{incompressible}(\\vec{\\lambda}_{incompressible})+\\kappa(J\\log{J} - J)``

Example Implementation: `StrainEnergyDensityFunction(LogarithmicCompressible(NeoHookean()), λ⃗, (ψ = (μ=μ, ), κ = κ))`
"""
struct LogarithmicCompressible <: AbstractHyperelasticModel
    incompressible::AbstractHyperelasticModel
end

function NonlinearContinua.StrainEnergyDensity(ψ::LogarithmicCompressible, λ⃗, p)
    StrainEnergyDensity(ψ.incompressible, λ⃗, p.ψ) + p.κ * (J * log(J) - J)
end

function NonlinearContinua.CauchyStressTensor(ψ::LogarithmicCompressible, λ⃗::Vector, p; adb=AD.ForwardDiffBackend())
    σ_dev = p.κ * (log(J(λ⃗)))
    σ = CauchyStressTensor(ψ.incompressible, λ⃗ ./ cbrt(prod(λ⃗)), p.ψ, adb=adb)
    return σ .+ σ_dev
end

function NonlinearContinua.CauchyStressTensor(ψ::LogarithmicCompressible, F::Matrix, p; adb=AD.ForwardDiffBackend())
    σ_dev = p.κ * (log(J(F)))
    σ = CauchyStressTensor(ψ.incompressible, F ./ cbrt(J(F)), p.ψ, adb=adb)
    return σ .+ σ_dev
end

function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::LogarithmicCompressible, λ⃗::Vector, p; adb=AD.ForwardDiffBackend())
    s_dev = p.κ * (log(J(λ))) ./ λ
    s = SecondPiolaKirchoffStressTensor(ψ.incompressible, λ⃗ ./ cbrt(J(λ⃗)), p.ψ, adb=adb)
    return s .+ s_dev
end

function NonlinearContinua.SecondPiolaKirchoffStressTensor(ψ::LogarithmicCompressible, F::Matrix, p; adb=AD.ForwardDiffBackend())
    s_dev = p.κ * (log(J(λ))) ./ λinv(F)
    s = SecondPiolaKirchoffStressTensor(ψ.incompressible, F ./ cbrt(J(F)), p.ψ, adb=adb)
    return s .+ s_dev
end
