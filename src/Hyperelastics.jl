"""
$(DocStringExtensions.README)
"""
module Hyperelastics

using DocStringExtensions

using PackageExtensionCompat
function __init__()
    @require_extensions
end

using Reexport
@reexport using ContinuumMechanicsBase
@reexport using InverseLangevinApproximations

using LossFunctions

using SpecialFunctions
using DataInterpolations
using QuadGK
using ComponentArrays, LabelledArrays, StructArrays
using LinearAlgebra, Statistics
# using SciMLBase

export HyperelasticUniaxialTest, HyperelasticBiaxialTest
export HyperelasticProblem
export predict
export parameters, parameter_bounds
export InvariantForm, PrincipalValueForm, DataDrivenForm

abstract type AbstractHyperelasticTest{T,S} <: ContinuumMechanicsBase.AbstractMaterialTest end
abstract type AbstractHyperelasticModel{T} <: ContinuumMechanicsBase.AbstractMaterialModel end

abstract type AbstractIncompressibleModel{T} <: AbstractHyperelasticModel{T} end
abstract type AbstractCompressibleModel{T} <: AbstractHyperelasticModel{T} end
abstract type AbstractDataDrivenHyperelasticModel{T} <: AbstractHyperelasticModel{T} end

# abstract type AbstractHyperelasticProblem <: AbstractOptimizationProblem end

struct InvariantForm end
struct PrincipalValueForm end
struct DataDrivenForm end

"""
$(SIGNATURES)
Creates an `OptimizationProblem` for use in [`Optimization.jl`](https://docs.sciml.ai/Optimization/stable/) to find the optimal parameters.

Fields:
- `ψ`: material model to use
- `test` or `tests`: A single or vector of hyperelastics tests to use when fitting the parameters
- `u₀`: Initial guess for parameters
- `ps`: Any additional parameters for calling predict
- `adb`: Select differentiation type from [`ADTypes.jl`](https://github.com/SciML/ADTypes.jl). The type is automatically applied to the type of AD applied to the Optimization Problem also.
- `loss`: Loss function from [`LossFunctions.jl`](https://github.com/JuliaML/LossFunctions.jl)
"""
function HyperelasticProblem end

# include("../ext/HyperelasticsOptimizationExt.jl")
include("invariants.jl")
include("material_tests.jl")
include("model_functions.jl")
include("stress_functions.jl")
include("datasets.jl")

include("general_form_incompressible_isotropic_models.jl")
include("isotropic_incompressible_models.jl")
include("isotropic_compressible_models.jl")

include("data_driven.jl")
include("average_chain_behavior.jl")

end
