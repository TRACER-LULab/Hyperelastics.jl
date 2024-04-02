# Hyperelastics
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://TRACER-LULab.github.io/Hyperelastics.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://TRACER-LULab.github.io/Hyperelastics.jl/dev)
[![codecov](https://codecov.io/gh/TRACER-LULab/Hyperelastics.jl/graph/badge.svg?token=EML9TQUEP9)](https://codecov.io/gh/TRACER-LULab/Hyperelastics.jl)


A hyperelastic model library and fitting toolkit developed by TRACER Lab at Liberty University. An extension is provided for [Optimization.jl](https://github.com/SciML/Optimization.jl) for model calibration based on experimental data.

## Installation:

To install `Hyperelastics.jl` in Julia >= v1.9, use the Julia package manager:

```julia
using Pkg
Pkg.add("Hyperelastics")
```

## Statement of Need:

The development of `Hyperelastics.jl` began as a study of the accuracy for a variety of material models for a set of experimental data. Often, researchers rely on custom implementations of material models and the data fitting process to find material parameters that match their experimental data. Hyperelastic models can well represent the nonlinear stress-deformation behavior of many biological tissues as well as engineering polymeric materials.

The SEDFs included in this package cover most (if not all) of the available analytical models from the literature to date, from constitutive to phenomelogical models. Furthermore, a selection of data-driven models are incldued as a starting point for the development of new methods.

`Hyperelastics.jl` is part of a spinoff Multi-Scale Material Modelling ($M^3$) Suite being developed by Vagus LLC (wwww.vagusllc.com), as a byproduct result of ongoing multi-functional material research being carried out in the Translational Robotics and Controls Engineering Research (TRACER) Lab at Liberty University. A pure Julia implementation allows for the use of automatic differentiation (AD) packages to calculate the partial derivatives of the SEDF. `Hyperelastics.jl` is designed to leverage multiple-dispatch to define a common set of functions for calculating the SED, Second Piola Kirchoff Stress Tensor, and the Cauchy Stress Tensor. The package provides a set of hyperelastic models and an interface to `Optimization.jl` for fitting model parameters. 

Currently, most commercial finite element codes only offer a limited number, often less than 10, of hyperelastic models which limits the extent to which researchers are able to accurately model a given material. The closest project to `Hyperelastics.jl` is the `matADi` project by Andreas Dutzler [@matAdi2023] which has AD support for 18 material models. 

## Citations:
All relevant citations are located in `CITATIONS.bib`

## Paper:
Markdown:
[![DOI](https://joss.theoj.org/papers/10.21105/joss.06314/status.svg)](https://doi.org/10.21105/joss.06314)

HTML:
<a style="border-width:0" href="https://doi.org/10.21105/joss.06314">
  <img src="https://joss.theoj.org/papers/10.21105/joss.06314/status.svg" alt="DOI badge" >
</a>

reStructuredText:
.. image:: https://joss.theoj.org/papers/10.21105/joss.06314/status.svg
   :target: https://doi.org/10.21105/joss.06314
