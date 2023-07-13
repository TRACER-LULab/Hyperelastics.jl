@testset "Material Tests" begin
    @test HyperelasticUniaxialTest(ones(3), ones(3), name = "test") isa
          HyperelasticUniaxialTest
    @test HyperelasticUniaxialTest(
        ones(3),
        ones(3),
        incompressible = false,
        name = "test",
    ) isa HyperelasticUniaxialTest

    @test HyperelasticUniaxialTest(ones(3), name = "test") isa HyperelasticUniaxialTest
    @test HyperelasticUniaxialTest(ones(3), incompressible = false, name = "test") isa
          HyperelasticUniaxialTest

    @test HyperelasticBiaxialTest(ones(3), ones(3), ones(3), ones(3), name = "test") isa
          HyperelasticBiaxialTest
    @test HyperelasticBiaxialTest(
        ones(3),
        ones(3),
        ones(3),
        ones(3),
        name = "test",
        incompressible = false,
    ) isa HyperelasticBiaxialTest

    @test HyperelasticBiaxialTest(ones(3), ones(3), name = "test") isa
          HyperelasticBiaxialTest
    @test HyperelasticBiaxialTest(
        ones(3),
        ones(3),
        name = "test",
        incompressible = false,
    ) isa HyperelasticBiaxialTest

    @test Treloar1944Uniaxial() isa HyperelasticUniaxialTest
    for i in [
        1.040,
        1.060,
        1.080,
        1.100,
        1.120,
        1.14,
        1.16,
        1.2,
        1.24,
        1.3,
        1.6,
        1.9,
        2.2,
        2.5,
        2.8,
        3.1,
        3.4,
        3.7,
    ]
        @test Kawabata1981(i) isa HyperelasticBiaxialTest
    end
end
