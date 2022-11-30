export GeneralMooneyRivlin, GeneralDarijaniNaghdabadi, GeneralBeda

"""
General Mooney Rivlin

Model:
```math
W = \\sum\\limits_{i,j = 0}^{N,M} C_{i,j}(I_1-3)^i(I_2-3)^j
```

Parameters:
- C: Matrix of coefficients for Cⱼᵢ coefficients

> Mooney M. A theory of large elastic deformation. Journal of applied physics. 1940 Sep;11(9):582-92.
"""
struct GeneralMooneyRivlin <: AbstractHyperelasticModel end

function NonlinearContinua.StrainEnergyDensity(ψ::GeneralMooneyRivlin, λ⃗::AbstractVector, (; C))
    I1 = I₁(λ⃗)
    I2 = I₂(λ⃗)
    I1_vec = map(Base.Fix1(^, I1), range(0, size(C, 2) - 1))
    I2_vec = map(Base.Fix1(^, I2), range(0, size(C, 1) - 1))
    W = I2_vec' * C * I1_vec
end

function NonlinearContinua.StrainEnergyDensity(ψ::GeneralMooneyRivlin, I⃗::AbstractVector, (; C), I::InvariantForm)
    @tullio W := C[j, i] * (I⃗[1] - 3)^(i - 1) * (I⃗[2] - 3)^(j - 1)
    return W
end

parameters(ψ::GeneralMooneyRivlin) = (:C,)

"""
General Darijani Naghdabadi

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
struct GeneralDarijaniNaghdabadi <: AbstractHyperelasticModel end

function NonlinearContinua.StrainEnergyDensity(ψ::GeneralDarijaniNaghdabadi, λ⃗::AbstractVector, (; A⃗, B⃗, m⃗, n⃗))
    @assert length(A⃗) == length(m⃗) "Length of A⃗ ≠ length of m⃗"
    @assert length(B⃗) == length(n⃗) "Length of B⃗ ≠ length of n⃗"
    sum(i -> sum(A⃗ .* (λ⃗[i] .^ m⃗ .- 1)) + sum(B⃗ .* (λ⃗[i] .^ (-1 .* n⃗) .- 1)), 1:3)
end

function parameters(ψ::GeneralDarijaniNaghdabadi)
    return (:A⃗, :B⃗, :m⃗, :n⃗)
end

"""
General Beda

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
struct GeneralBeda <: AbstractHyperelasticModel end

function NonlinearContinua.StrainEnergyDensity(ψ::GeneralBeda, λ⃗::AbstractVector, (; C, K, α, β))
    @assert length(C) == length(α) "Vector C and Vector α are not the same length"
    @assert length(K) == length(β) "Vector K and Vector β are not the same length"
    @tullio W1 := C[i] / α[i] * (I₁(λ⃗) - 3)^α[i]# |> sum
    @tullio W2 := K[i] / β[i] * (I₂(λ⃗) - 3)^β[i]# |> sum
    W = W1 + W2
    return W
end

function NonlinearContinua.StrainEnergyDensity(ψ::GeneralBeda, I⃗::AbstractVector, (; C, K, α, β), I::InvariantForm)
    @assert length(C) == length(α) "Vector C and Vector α are not the same length"
    @assert length(K) == length(β) "Vector K and Vector β are not the same length"
    W1 = C ./ α .* (I⃗[1] - 3) .^ α |> sum
    W2 = K ./ β .* (I⃗[2] - 3) .^ β |> sum
    return W1 + W2
end

parameters(ψ::GeneralBeda) = (:C, :K, :α, :β)
