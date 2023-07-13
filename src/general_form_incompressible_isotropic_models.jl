export GeneralMooneyRivlin, GeneralDarijaniNaghdabadi, GeneralBeda

"""
`GeneralMooneyRivlin`

Model:

```math
W = \\sum\\limits_{i,j = 0}^{N,M} C_{i,j}(I_1-3)^i(I_2-3)^j
```

Parameters:
- C: Matrix of coefficients for Cⱼᵢ coefficients

> Mooney M. A theory of large elastic deformation. Journal of applied physics. 1940 Sep;11(9):582-92.
"""
struct GeneralMooneyRivlin{T} <: AbstractIncompressibleModel{T}
    f::Function
    function GeneralMooneyRivlin(
        I::Union{InvariantForm,PrincipalValueForm} = PrincipalValueForm(),
    )
        function f(x, (; I1, I2, C))
            return C[x] * I1^(x[2] - 1) * I2^(x[1] - 1)
        end
        new{typeof(I)}(f)
    end
end

function NonlinearContinua.StrainEnergyDensity(
    ψ::GeneralMooneyRivlin{I},
    λ⃗::Vector{T},
    (; C⃗),
) where {I<:PrincipalValueForm,T}
    I1 = I₁(λ⃗)
    I2 = I₂(λ⃗)
    W = sum(Base.Fix2(ψ.f, (I1 = I1, I2 = I2, C = C⃗)), CartesianIndices(C⃗))
    # W = [C[j+1, i+1] * I1^i * I2^j for i in 0:size(C, 2)-1 for j in 0:size(C, 1)-1] |> sum

    return W
end

function NonlinearContinua.StrainEnergyDensity(
    ψ::GeneralMooneyRivlin{I},
    I⃗::Vector{T},
    (; C⃗),
) where {I<:InvariantForm,T}
    W = sum(Base.Fix2(ψ.f, (I1 = I⃗[1], I2 = I⃗[2], C = C⃗)), CartesianIndices(C⃗))
    return W
end

parameters(::GeneralMooneyRivlin) = (:C⃗,)

"""
`GeneralDarijaniNaghdabadi`

Model:

```math
W = \\sum\\limits_{i = 1}^{3}\\sum\\limits_{j=0}^{N} A_j (\\lambda_i^{m_j}-1) + B_j(\\lambda_i^{-n_j}-1)
```

Parameters:
- A⃗
- B⃗
- m⃗
- n⃗

> Bahreman M, Darijani H. New polynomial strain energy function; application to rubbery circular cylinders under finite extension and torsion. Journal of Applied Polymer Science. 2015 Apr 5;132(13).
"""
struct GeneralDarijaniNaghdabadi{T} <: AbstractIncompressibleModel{T}
    GeneralDarijaniNaghdabadi(::T = PrincipalValueForm()) where {T<:PrincipalValueForm} =
        new{T}()
end

function NonlinearContinua.StrainEnergyDensity(
    ::GeneralDarijaniNaghdabadi,
    λ⃗::Vector{T},
    (; A⃗, B⃗, m⃗, n⃗),
) where {T}

    w1 = sum(A⃗ .* (λ⃗[1] .^ m⃗ .- 1)) + sum(B⃗ .* (λ⃗[1] .^ (-1 .* n⃗) .- 1))
    w2 = sum(A⃗ .* (λ⃗[2] .^ m⃗ .- 1)) + sum(B⃗ .* (λ⃗[2] .^ (-1 .* n⃗) .- 1))
    w3 = sum(A⃗ .* (λ⃗[3] .^ m⃗ .- 1)) + sum(B⃗ .* (λ⃗[3] .^ (-1 .* n⃗) .- 1))

    return w1 + w2 + w3
end

function parameters(::GeneralDarijaniNaghdabadi)
    return (:A⃗, :B⃗, :m⃗, :n⃗)
end

"""
`GeneralBeda`

Model:

```math
W = \\sum\\limits_{i = 1}^{N}\\frac{C_i}{\\alpha_i}(I_1-3)^{\\alpha_i} + \\sum\\limits_{j=1}^{M}\\frac{K_j}{\\beta_j}(I_2-3)^{\\beta_j}
```

Parameters:
- C
- K
- α
- β

> Beda T. Reconciling the fundamental phenomenological expression of the strain energy of rubber with established experimental facts. Journal of Polymer Science Part B: Polymer Physics. 2005 Jan 15;43(2):125-34.
"""
struct GeneralBeda{T} <: AbstractIncompressibleModel{T}
    GeneralBeda(I::Union{InvariantForm,PrincipalValueForm} = PrincipalValueForm()) =
        new{typeof(I)}()
end

function NonlinearContinua.StrainEnergyDensity(
    ::GeneralBeda{I},
    λ⃗::Vector{T},
    (; C⃗, K⃗, α⃗, β⃗),
) where {I<:PrincipalValueForm,T}
    @assert length(C⃗) == length(α⃗) "Vector C and Vector α are not the same length"
    @assert length(K⃗) == length(β⃗) "Vector K and Vector β are not the same length"
    I1 = I₁(λ⃗)
    I2 = I₂(λ⃗)
    W1 = @. C⃗ / α⃗ * (I1 - 3)^α⃗ # |> sum
    W2 = @. K⃗ / β⃗ * (I2 - 3)^β⃗ # |> sum
    W = sum(W1) + sum(W2)
    return W
end

function NonlinearContinua.StrainEnergyDensity(
    ::GeneralBeda{I},
    I⃗::Vector{T},
    (; C⃗, K⃗, α⃗, β⃗),
) where {I<:InvariantForm,T}
    @assert length(C⃗) == length(α⃗) "Vector C and Vector α are not the same length"
    @assert length(K⃗) == length(β⃗) "Vector K and Vector β are not the same length"
    W1 = @. C⃗ / α⃗ * (I⃗[1] - 3)^α⃗
    W2 = @. K⃗ / β⃗ * (I⃗[2] - 3)^β⃗
    return sum(W1) + sum(W2)
end

parameters(::GeneralBeda) = (:C⃗, :K⃗, :α⃗, :β⃗)
