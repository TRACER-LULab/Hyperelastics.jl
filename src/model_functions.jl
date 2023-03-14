"""
`parameters(ψ::AbstractHyperelasticModel)`

Returns a tuple of the parameters required for the model

Fields
- `ψ`: Hyperelastics model
"""
function parameters(ψ::AbstractHyperelasticModel)
    @error "$(typeof(ψ)) does not have a parameters function implemented"
end

"""
`parameter_bounds(ψ::AbstractHyperelasticModel, test::AbstractHyperelasticTest)`
`parameter_bounds(ψ::AbstractHyperelasticModel, tests::Vector{AbstractHyperelasticTest})`

Returns a tuple of the parameter bounds provided the experimental data and model

Fields
- `ψ`: Hyperelastic model
- `test` or `tests`: The test or vector of tests to use in finding the parameter bounds.
"""
function parameter_bounds(ψ::AbstractHyperelasticModel, test::AbstractHyperelasticTest)
    lb = nothing
    ub = nothing
    return (lb=lb, ub=ub)
end

function parameter_bounds(ψ::AbstractHyperelasticModel, tests::Vector{<:AbstractHyperelasticTest})
    bounds = map(Base.Fix1(parameter_bounds, ψ), tests)
    lbs = getfield.(bounds, :lb)
    ubs = getfield.(bounds, :ub)
    if !(eltype(lbs) <: Nothing)
        lb_ps = fieldnames(eltype(lbs))
        lb = map(p -> p .=> maximum(getfield.(lbs, p)), lb_ps) |> NamedTuple
    else
        lb = nothing
    end

    if !(eltype(ubs) <: Nothing)
        ub_ps = fieldnames(eltype(ubs))
        ub = map(p -> p .=> minimum(getfield.(ubs, p)), ub_ps) |> NamedTuple
    else
        ub = nothing
    end
    return (lb=lb, ub=ub)
end

"""
Returns a constraint equation for models were parameters bounds are interdependent
"""
function constraints(ψ::AbstractHyperelasticModel, data::AbstractHyperelasticTest)
    return nothing
end
