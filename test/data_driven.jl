@testset "Data-Driven models" begin
    λ⃗ = [1.1, inv(sqrt(1.1)), inv(sqrt(1.1))]
    F = diagm(λ⃗)

    data_uniaxial = Treloar1944Uniaxial()

    λ₁ = [1.04, 3.1]
    data_biaxial = Kawabata1981.(λ₁)

    λ⃗ = [1.1, inv(sqrt(1.1)), inv(sqrt(1.1))]
    F = diagm(λ⃗)

    for model in [SussmanBathe, DataDrivenAverageChainBehavior]#, MacroMicroMacro]
        # if model == MacroMicroMacro
        #     n1 = 100
        #     p₀ = range(0.0, 0.7, length=n1) |> collect
        #     λ_max = maximum(maximum.(map(x -> maximum.(x.data.λ), tests)))
        #     λs = collect(range(0.001, 6.0, length=n1))
        #     PChain(u) = BSplineApprox(u, λs, 3, 12, :Uniform, :Uniform)
        #     ψ = model(data_biaxial, PChain, p₀)
        # else
            ψ = model(data_uniaxial)
        # end
        for deformation in [λ⃗, F]
            W = StrainEnergyDensity(ψ, deformation, [])
            S = SecondPiolaKirchoffStressTensor(ψ, deformation, nothing)
            σ = CauchyStressTensor(ψ, deformation, nothing)
            @test sum(isnan.(W)) == 0
            @test sum(isinf.(W)) == 0
            @test sum(isnan.(S)) == 0
            @test sum(isinf.(S)) == 0
            @test sum(isnan.(σ)) == 0
            @test sum(isinf.(σ)) == 0
        end
    end
end
