"""
UniaxialHyperelasticData(s₁, λ₁)

Create a uniaxial hyperelastic data object from arrays of test data. The function returns a HyperelasticData object with the stresses and principal stretches. Currently, this assumes the material is incompressible.
"""
struct UniaxialHyperelasticData{T} <: AbstractHyperelasticData
    s⃗::T
    λ⃗::T
    function UniaxialHyperelasticData(s₁, λ₁)
        @assert length(s₁) == length(λ₁) "Length of stress must be equal to the length of the stretches"
        f(x) = [x, 1/√x, 1/√x]
        λ⃗ = map(f, λ₁)
        g(x) = [x]
        s⃗ = map(g, s₁)
        return new{promote_type(typeof(s⃗), typeof(λ⃗))}(s⃗, λ⃗)
    end
end

struct BiaxialHyperelasticData{T} <: AbstractHyperelasticData
    s⃗::T
    λ⃗::T
end

"""
BiaxialHyperelasticData(s₁, s₂, λ₁, λ₂)

Create a biaxial hyperelastic data object from arrays of test data. The function returns a HyperelasticData object with the stresses and principal stretches. Currently, this assumes the material is incompressible.
"""
function BiaxialHyperelasticData(s₁, s₂, λ₁, λ₂)
    f(x) = [x[1], x[2], 1/x[1]/x[2]]
    λ⃗ = map(f, zip(λ₁, λ₂))
    g(x) = [x[1], x[2]]
    s⃗ = map(g, zip(s₁, s₂))
    return BiaxialHyperelasticData{promote_type(typeof(s⃗), typeof(λ⃗))}(s⃗, λ⃗)
end
