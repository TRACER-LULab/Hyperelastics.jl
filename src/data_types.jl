"""
UniaxialHyperelasticData(s₁, λ₁)

Create a uniaxial hyperelastic data object from arrays of test data. The function returns a HyperelasticData object with the stresses and principal stretches. Currently, this assumes the material is incompressible.
"""
struct UniaxialHyperelasticData <: AbstractHyperelasticData
    s⃗
    λ⃗
    UniaxialHyperelasticData(s⃗, λ⃗) =
        let
            s⃗ = zip(s⃗)
            λ⃗ = zip(λ⃗, (λ⃗) .^ (-0.5), (λ⃗) .^ (-0.5))
            new(s⃗, λ⃗)
        end
end

struct BiaxialHyperelasticData <: AbstractHyperelasticData
    s⃗
    λ⃗
end

"""
BiaxialHyperelasticData(s₁, s₂, λ₁, λ₂)

Create a biaxial hyperelastic data object from arrays of test data. The function returns a HyperelasticData object with the stresses and principal stretches. Currently, this assumes the material is incompressible.
"""
function BiaxialHyperelasticData(s₁, s₂, λ₁, λ₂)
    s⃗ = zip(s₁, s₂)
    λ⃗ = zip(λ₁, λ₂, (λ₁ .* λ₂) .^ (-1))
    return BiaxialHyperelasticData(s⃗, λ⃗)
end
