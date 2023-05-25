struct HyperelasticDataEntry{T,S}
    λ::Vector{T}
    s::Vector{S}
end

"""
`HyperelasticUniaxialTest(λ₁, s₁; name, incompressible=true)`
`HyperelasticUniaxialTest(λ₁; name, incompressible=true)`

Creates an object storing results from a uniaxial test of a hyperelatic  material.

Fields:
- `λ₁`: Vector of uniaxial stretches
- `s₁`: Vector of experiemntal stresses (optional)
- `name`: string for the name of the test
- `incompressible`: `true` if the material can be assumed to be incompressible.
"""
struct HyperelasticUniaxialTest{T,S} <: AbstractHyperelasticTest{T,S}
    data::StructVector
    name::String
    function HyperelasticUniaxialTest(λ₁, s₁; name, incompressible=true)
        @assert length(λ₁) == length(s₁) "Inputs must be the same length"
        if incompressible
            λ₂ = λ₃ = @. sqrt(1 / λ₁)
        else
            λ₂ = λ₃ = Vector{eltype(λ₁)}(undef, length(λ₁))
        end
        λ = collect.(zip(λ₁, λ₂, λ₃))
        s = collect.(zip(s₁))
        data = StructArray{HyperelasticDataEntry}((λ, s))
        new{eltype(eltype(λ)), eltype(eltype(s))}(data, name)
    end
    function HyperelasticUniaxialTest(λ₁; name, incompressible=true)
        if incompressible
            λ₂ = λ₃ = @. sqrt(1 / λ₁)
        else
            λ₂ = λ₃ = Vector{eltype(λ₁)}(undef, length(λ₁))
        end
        λ = collect.(zip(λ₁, λ₂, λ₃))
        s = collect.(zip(Vector{eltype(λ₁)}(undef, length(λ₁))))
        data = StructArray{HyperelasticDataEntry}((λ, s))
        new{eltype(eltype(λ)), eltype(eltype(s))}(data, name)
    end
end

function Base.show(io::IO, test::HyperelasticUniaxialTest)
    print(io, Term.RenderableText("Uniaxial Test: {bold} $(test.name)"))
    print(io,
        Term.Table(
            hcat(getindex.(test.data.λ, 1), getindex.(test.data.s, 1)),
            header=["λ₁", "s₁"],
            box=:ROUNDED,
        )
    )
end

"""
`HyperelasticBiaxialTest(λ₁, λ₂, s₁, s₂; name, incompressible=true)`
`HyperelasticBiaxialTest(λ₁, λ₂; name, incompressible=true)`

Creates an object storing results from a biaxial test of a hyperelatic material.

Fields:
- `λ₁`: Vector of 1-direction stretches
- `λ₂`: Vector of 2-direction stretchs
- `s₁`: Vector of experiemntal stresses (optional)
- `s₂`: Vector of experiemntal stresses (optional)
- `name`: string for the name of the test
- `incompressible`: `true` if the material can be assumed to be incompressible.
"""
struct HyperelasticBiaxialTest{T,S} <: AbstractHyperelasticTest{T,S}
    data::StructVector
    name::String
    function HyperelasticBiaxialTest(λ₁, λ₂, s₁, s₂; name, incompressible=true)
        @assert length(λ₁) == length(λ₂) == length(s₁) == length(s₂) "Inputs must be the same length"
        if incompressible
            λ₃ = @. 1 / λ₁ / λ₂
        else
            λ₃ = Vector{eltype(λ₁)}(undef, length(λ₁))
        end
        λ = collect.(zip(λ₁, λ₂, λ₃))
        s = collect.(zip(s₁, s₂))
        data = StructArray{HyperelasticDataEntry}((λ, s))
        new{promote(eltype(eltype(λ₁)), eltype(eltype(λ₂))),promote(eltype(eltype(s₁)), eltype(eltype(s₂)))}(data, name)
    end
    function HyperelasticBiaxialTest(λ₁, λ₂; name, incompressible=true)
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
        new{promote(eltype(eltype(λ₁)), eltype(eltype(λ₂))),promote(eltype(eltype(s₁)), eltype(eltype(s₂)))}(data, name)
    end
end

function Base.show(io::IO, test::HyperelasticBiaxialTest)
    print(io, Term.RenderableText("Biaxial Test: {bold} $(test.name)"))
    print(io,
        Term.Table(
            hcat(getindex.(test.data.λ, 1), getindex.(test.data.λ, 2), getindex.(test.data.s, 1), getindex.(test.data.s, 2)),
            header=["λ₁", "λ₂", "s₁", "s₂"],
            box=:ROUNDED,
        )
    )
end

function NonlinearContinua.predict(ψ::AbstractHyperelasticModel, test::HyperelasticUniaxialTest{T,S}, p, adtype; kwargs...) where {T,S}
    f(λ) = SecondPiolaKirchoffStressTensor(ψ, λ, p, adtype, return_type=S, kwargs...)
    λ = test.data.λ
    s = map(f, λ)
    s₁ = getindex.(s, 1)
    s₃ = getindex.(s, 3)
    λ₁ = getindex.(λ, 1)
    λ₃ = getindex.(λ, 3)
    Δs₁₃ = @. s₁ - s₃ * λ₃ / λ₁
    pred = HyperelasticUniaxialTest(λ₁, Δs₁₃, name=test.name)
    return pred
end

function NonlinearContinua.predict(ψ::AbstractHyperelasticModel, test::HyperelasticBiaxialTest, p, adtype; kwargs...)
    f(λ) = SecondPiolaKirchoffStressTensor(ψ, λ, p, adtype, kwargs...)
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
    pred = HyperelasticBiaxialTest(λ₁, λ₂, Δs₁₃, Δs₂₃, name=test.name)
    return pred
end
