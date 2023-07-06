@testset "Material Tests" begin
    @test HyperelasticUniaxialTest(ones(3), ones(3), name="test") isa HyperelasticUniaxialTest
    @test HyperelasticUniaxialTest(ones(3), ones(3), incompressible=false,name="test") isa HyperelasticUniaxialTest

    @test HyperelasticUniaxialTest(ones(3), name="test") isa HyperelasticUniaxialTest
    @test HyperelasticUniaxialTest(ones(3), incompressible=false, name="test") isa HyperelasticUniaxialTest

    @test HyperelasticBiaxialTest(ones(3), ones(3), ones(3), ones(3), name="test") isa HyperelasticBiaxialTest
    @test HyperelasticBiaxialTest(ones(3), ones(3), ones(3), ones(3), name="test", incompressible=false) isa HyperelasticBiaxialTest

    @test HyperelasticBiaxialTest(ones(3), ones(3), name="test", ) isa HyperelasticBiaxialTest
    @test HyperelasticBiaxialTest(ones(3), ones(3), name="test", incompressible=false) isa HyperelasticBiaxialTest
end
