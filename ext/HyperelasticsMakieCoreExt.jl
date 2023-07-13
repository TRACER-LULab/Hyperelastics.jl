module HyperelasticsMakieCoreExt

using MakieCore
using Hyperelastics

function MakieCore.convert_arguments(
    p::MakieCore.PointBased,
    test::HyperelasticUniaxialTest,
)
    MakieCore.convert_arguments(p, getindex.(test.data.λ, 1), getindex.(test.data.s, 1))
end

end
