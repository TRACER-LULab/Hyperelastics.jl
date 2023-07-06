using LinearAlgebra
using Hyperelastics
include("../../NonlinearContinua/src/NonlinearContinua.jl")
using Test
using ForwardDiff, FiniteDiff, Zygote, Enzyme
using InteractiveUtils

@time include("helpers.jl")
@time include("material_tests.jl")
@time include("compressible_models.jl")
@time include("data_driven.jl")
