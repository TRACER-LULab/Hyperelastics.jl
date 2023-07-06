@testset "Helper Functions" begin
    λ⃗ = [1.1, inv(sqrt(1.1)), inv(sqrt(1.1))]
    F = diagm(λ⃗)

    @test I₁(λ⃗) isa Float64
    @test I₁(F) isa Float64
    @test I₂(λ⃗) isa Float64
    @test I₂(F) isa Float64
    @test I₃(λ⃗) isa Float64
    @test I₃(F) isa Float64
    @test J(λ⃗) isa Float64
    @test J(F) isa Float64

end
