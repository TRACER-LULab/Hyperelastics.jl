```@meta
CurrentModule = Hyperelastics
```

# Hyperelastics

`Hyperelastics.jl` relies on `ContinuumMechanicsBase.jl` for the definitions of key continuum terms. Please cite the package if you use it in your work.

## Installation 

To install `Hyperelastics.jl`, use the Julia package manager:

```julia
using Pkg
Pkg.add("Hyperelastics")
```

## Statement of Need

The development of `Hyperelastics.jl` began as a study of the accuracy for a variety of material models for a set of experimental data. Often, researchers rely on custom implementations of material models and the data fitting process to find material parameters that match their experimental data. Hyperelastic models can well represent the nonlinear stress-deformation behavior of many biological tissues as well as engineering polymeric materials.

The SEDFs included in this package cover most (if not all) of the available analytical models from the literature to date, from constitutive to phenomelogical models. Furthermore, a selection of data-driven models are incldued as a starting point for the development of new methods.

`Hyperelastics.jl` is part of a spinoff Multi-Scale Material Modelling ($M^3$) Suite being developed by Vagus LLC (www.vagusllc.com), as a byproduct result of ongoing multi-functional material research being carried out in the Translational Robotics and Controls Engineering Research (TRACER) Lab at Liberty University. A pure Julia implementation allows for the use of automatic differentiation (AD) packages to calculate the partial derivatives of the SEDF. `Hyperelastics.jl` is designed to leverage multiple-dispatch to define a common set of functions for calculating the SED, Second Piola Kirchoff Stress Tensor, and the Cauchy Stress Tensor. The package provides a set of hyperelastic models and an interface to `Optimization.jl` for fitting model parameters. 

Currently, most commercial finite element codes only offer a limited number, often less than 10, of hyperelastic models which limits the extent to which researchers are able to accurately model a given material. The closest project to `Hyperelastics.jl` is the `matADi` project by Andreas Dutzler [@matAdi2023] which has AD support for 18 material models. 

## Community Guidelines

To sustainably develop the package, we will use the established practices from the SciML community for guidlines:
- For any issues with or contributions, please open an [issue in GitHub](https://github.com/TRACER-LULab/Hyperelastics.jl/issues) or a [PR](https://github.com/TRACER-LULab/Hyperelastics.jl/pulls) as appropriate.

- Please refer to the
    [SciML ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://github.com/SciML/ColPrac/blob/master/README.md)
    for guidance on PRs, issues, and other matters relating to contributing.

- See the [SciML Style Guide](https://github.com/SciML/SciMLStyle) for common coding practices and other style decisions.
