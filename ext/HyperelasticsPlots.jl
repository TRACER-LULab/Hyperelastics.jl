module HyperelasticsPlots

import Plots: Plots, @ext_imp_use, @recipe, @series
@ext_imp_use :import Hyperelastics HyperelasticUniaxialTest HyperelasticBiaxialTest

@recipe function f(test::HyperelasticUniaxialTest)
    @series begin
        seriestype := :scatter
        getindex.(test.data.λ, 1),getindex.(test.data.s, 1)
    end
    # getindex.(test.data.λ, 1),getindex.(test.data.s, 1)
end

@recipe function f(test::HyperelasticBiaxialTest)
    xlabel --> "Stretch"
    ylabel --> "Stress"
    @series begin
        seriestype := :scatter
        label := "Stretch 1 - Stress 1"
        getindex.(test.data.λ, 1),getindex.(test.data.s, 1)
    end

    @series begin
        seriestype := :scatter
        label := "Stretch 1 - Stress 2"
        getindex.(test.data.λ, 1),getindex.(test.data.s, 2)
    end

    @series begin
        seriestype := :scatter
        label := "Stretch 2 - Stress 1"
        getindex.(test.data.λ, 2),getindex.(test.data.s, 1)
    end

    @series begin
        seriestype := :scatter
        label := "Stretch 2 - Stress 2"
        getindex.(test.data.λ, 2),getindex.(test.data.s, 2)
    end
end

end
