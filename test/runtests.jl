using LinearAlgebra

using Hyperelastics
using Test
using ForwardDiff, FiniteDiff, Zygote, Enzyme
using InteractiveUtils
using ADTypes

@time include("helpers.jl")
@time include("material_tests.jl")
@time include("models.jl")
@time include("data_driven.jl")
@time include("model_fitting.jl")
