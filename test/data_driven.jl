@testset "Data-Driven models" begin
    λ⃗ = [1.1, inv(sqrt(1.1)), inv(sqrt(1.1))]
    F = diagm(λ⃗)

    data_uniaxial = Treloar1944Uniaxial()

    λ₁ = [1.04, 3.1]
    data_biaxial = Kawabata1981.(λ₁)

    λ⃗ = [1.1, inv(sqrt(1.1)), inv(sqrt(1.1))]
    F = diagm(λ⃗)

    for model in [SussmanBathe, DataDrivenAverageChainBehavior]
        @show model
        ψ = model(data_uniaxial)
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
