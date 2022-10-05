using Hyperelastics
using SnoopCompile, ProfileView

tinf = @snoopi_deep StrainEnergyDensityFunction(Gent(), [1.0, 1.0, 1.0], (μ=10.0,Jₘ=10.0))

ProfileView.view(flamegraph(tinf))
