struct HyperelasticDataEntry{T,S}
    λ::Vector{T}
    s::Vector{S}
end

struct HyperelasticUniaxialTest{T,S} <: AbstractHyperelasticTest{T,S}
    data::StructVector
    name::String
    """
    $(SIGNATURES)

    Creates an object storing results from a uniaxial test of a hyperelatic  material.

    # Arguments:
    - `λ₁`: Vector of uniaxial stretches
    - `s₁`: Vector of experiemntal stresses (optional)
    - `name`: string for the name of the test
    - `incompressible`: `true` if the material can be assumed to be incompressible.
    """
    function HyperelasticUniaxialTest(λ₁, s₁; name, incompressible = true)
        @assert length(λ₁) == length(s₁) "Inputs must be the same length"
        if incompressible
            λ₂ = λ₃ = @. sqrt(1 / λ₁)
        else
            λ₂ = λ₃ = Vector{eltype(λ₁)}(undef, length(λ₁))
        end
        λ = collect.(zip(λ₁, λ₂, λ₃))
        s = collect.(zip(s₁))
        data = StructArray{HyperelasticDataEntry}((λ, s))
        new{eltype(eltype(λ)),eltype(eltype(s))}(data, name)
    end
    function HyperelasticUniaxialTest(λ₁; name, incompressible = true)
        if incompressible
            λ₂ = λ₃ = @. sqrt(1 / λ₁)
        else
            λ₂ = λ₃ = Vector{eltype(λ₁)}(undef, length(λ₁))
        end
        λ = collect.(zip(λ₁, λ₂, λ₃))
        s = collect.(zip(Vector{eltype(λ₁)}(undef, length(λ₁))))
        data = StructArray{HyperelasticDataEntry}((λ, s))
        new{eltype(eltype(λ)),eltype(eltype(s))}(data, name)
    end
end

struct HyperelasticBiaxialTest{T,S} <: AbstractHyperelasticTest{T,S}
    data::StructVector
    name::String
    """
    $(SIGNATURES)

    Creates an object storing results from a biaxial test of a hyperelatic material.

    # Arguments:
    - `λ₁`: Vector of 1-direction stretches
    - `λ₂`: Vector of 2-direction stretchs
    - `s₁`: Vector of experiemntal stresses (optional)
    - `s₂`: Vector of experiemntal stresses (optional)
    - `name`: string for the name of the test
    - `incompressible`: `true` if the material can be assumed to be incompressible.
    """
    function HyperelasticBiaxialTest(λ₁, λ₂, s₁, s₂; name, incompressible = true)
        @assert length(λ₁) == length(λ₂) == length(s₁) == length(s₂) "Inputs must be the same length"
        if incompressible
            λ₃ = @. 1 / λ₁ / λ₂
        else
            λ₃ = Vector{eltype(λ₁)}(undef, length(λ₁))
        end
        λ = collect.(zip(λ₁, λ₂, λ₃))
        s = collect.(zip(s₁, s₂))
        data = StructArray{HyperelasticDataEntry}((λ, s))
        T = promote_type(eltype(eltype(λ₁)), eltype(eltype(λ₂)))
        S = promote_type(eltype(eltype(s₁)), eltype(eltype(s₂)))
        new{T,S}(data, name)
    end
    function HyperelasticBiaxialTest(λ₁, λ₂; name, incompressible = true)
        @assert length(λ₁) == length(λ₂) "Inputs must be the same length"
        if incompressible
            λ₃ = @. 1 / λ₁ / λ₂
        else
            λ₃ = Vector{eltype(λ₁)}(undef, length(λ₁))
        end
        s₁ = s₂ = Vector{eltype(λ₁)}(undef, length(λ₁))
        λ = collect.(zip(λ₁, λ₂, λ₃))
        s = collect.(zip(s₁, s₂))
        data = StructArray{HyperelasticDataEntry}((λ, s))
        new{
            promote_type(eltype(eltype(λ₁)), eltype(eltype(λ₂))),
            promote_type(eltype(eltype(s₁)), eltype(eltype(s₂))),
        }(
            data,
            name,
        )
    end
end

## Predict overloads
function ContinuumMechanicsBase.predict(
    ψ::AbstractHyperelasticModel{A},
    test::HyperelasticUniaxialTest{B,C},
    p;
    kwargs...,
) where {A,B,C}
    f(λ) = SecondPiolaKirchoffStressTensor(ψ, λ, p; kwargs...)
    λ = test.data.λ
    s = map(f, λ)
    s₁ = getindex.(s, 1)
    s₃ = getindex.(s, 3)
    λ₁ = getindex.(λ, 1)
    λ₃ = getindex.(λ, 3)
    Δs₁₃ = @. s₁ - s₃ * λ₃ / λ₁
    pred = HyperelasticUniaxialTest(λ₁, Δs₁₃, name = test.name)
    return pred
end

function ContinuumMechanicsBase.predict(
    ψ::AbstractHyperelasticModel{A},
    test::HyperelasticBiaxialTest{B,C},
    p;
    kwargs...,
) where {A,B,C}
    f(λ) = SecondPiolaKirchoffStressTensor(ψ, λ, p; kwargs...)
    λ = test.data.λ
    s = map(f, λ)
    s₁ = getindex.(s, 1)
    s₂ = getindex.(s, 2)
    s₃ = getindex.(s, 3)
    λ₁ = getindex.(λ, 1)
    λ₂ = getindex.(λ, 2)
    λ₃ = getindex.(λ, 3)
    Δs₂₃ = @. s₂ - s₃ * λ₃ / λ₁
    Δs₁₃ = @. s₁ - s₃ * λ₃ / λ₁
    pred = HyperelasticBiaxialTest(λ₁, λ₂, Δs₁₃, Δs₂₃, name = test.name)
    return pred
end
