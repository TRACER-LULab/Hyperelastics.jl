module Hyperelastics

using Reexport
@reexport using NonlinearContinua
using InverseLangevinApproximations
using LossFunctions, Optimization
using AbstractDifferentiation, ForwardDiff
using Tullio
using SpecialFunctions
using DataInterpolations
using QuadGK
using ComponentArrays, LabelledArrays, StructArrays
using LinearAlgebra, Statistics
using Term

export HyperelasticUniaxialTest, HyperelasticBiaxialTest, HyperelasticProblem
export predict
export citation, parameters, parameter_bounds

abstract type AbstractHyperelasticTest <: NonlinearContinua.AbstractMaterialTest end
abstract type AbstractHyperelasticModel <: NonlinearContinua.AbstractMaterialModel end
abstract type AbstractDataDrivenHyperelasticModel <: AbstractHyperelasticModel end
abstract type AbstractHyperelasticProblem end

struct InvariantForm end

include("data_types.jl")
include("invariants.jl")
include("model_functions.jl")
include("datasets.jl")
include("isotropic_incompressible_models.jl")
include("isotropic_compressible_models.jl")
include("data_driven.jl")
include("macro_micro_macro_model.jl")
include("average_chain_behavior.jl")
# include("optimization_interface.jl")
include("material_tests.jl")
# using SnoopPrecompile
# @precompile_setup begin
#     λ₁ = 2.0
#     data = UniaxialHyperelasticData([0.0], [λ₁])
#     ψ⃗ = filter(x -> typeof(x) <: DataType, Base.Fix1(getfield, Hyperelastics).(names(Hyperelastics)))
#     ψ⃗ = filter(x -> !(x <: AbstractDataDrivenHyperelasticModel) && (x <: AbstractHyperelasticModel), ψ⃗)

#     bounds = map(ψ -> parameter_bounds(ψ(), data), ψ⃗)

#     p⃗ = map(ψ -> NamedTuple(zip(parameters(ψ()), 2.0 * ones(length(parameters(ψ()))))), ψ⃗)

#     for i in eachindex(p⃗)
#         if !(isnothing(bounds[i].lb))
#             for (key, value) in pairs(bounds[i].lb)
#                 if p⃗[i][key] < value
#                     dict_form = Dict(pairs(p⃗[i]))
#                     dict_form[key] = value * 2
#                     p⃗[i] = NamedTuple(dict_form)
#                 end
#             end
#         end
#         if !(isnothing(bounds[i].ub))
#             for (key, value) in pairs(bounds[i].ub)
#                 if p⃗[i][key] > value
#                     dict_form = Dict(pairs(p⃗[i]))
#                     dict_form[key] = value / 1.5
#                     p⃗[i] = NamedTuple(dict_form)
#                 end
#             end
#         end
#         if ψ⃗[i] == Hyperelastics.GeneralMooneyRivlin
#             p⃗[i] = (C=ones(3, 3),)
#         end
#         if ψ⃗[i] == Hyperelastics.KhiemItskov
#             p⃗[i] = (μcκ=1.0, n=1.67, q=1.0, μt=1.0)
#         end
#         if ψ⃗[i] == Hyperelastics.Shariff
#             p⃗[i] = (E=1.0, α⃗=ones(5))
#         end
#         if ψ⃗[i] == Hyperelastics.VanDerWaals
#             p⃗[i] = (μ=10.0, λm=2 * sqrt(3), β=0.9, α=10.0)
#         end
#         if ψ⃗[i] == Hyperelastics.EdwardVilgis
#             p⃗[i] = (Ns=20, Nc=20, α=0.01, η=0.05)
#         end
#         @precompile_all_calls begin
#             p⃗_LabelledArrays = map(LVector, p⃗)
#             p⃗_ComponentArrays = map(ComponentVector, p⃗)
#             for (ψ, p) in zip(ψ⃗, p⃗)
#                 StrainEnergyDensity(ψ(), data.λ⃗[1], p)
#             end
#             for (ψ, p) in zip(ψ⃗, p⃗_LabelledArrays)
#                 StrainEnergyDensity(ψ(), data.λ⃗[1], p)
#             end
#             for (ψ, p) in zip(ψ⃗, p⃗_ComponentArrays)
#                 StrainEnergyDensity(ψ(), data.λ⃗[1], p)
#             end
#         end
#     end
# end

end
