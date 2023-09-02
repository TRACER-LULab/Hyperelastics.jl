---
title: 'Hyperelastics.jl: A Julia package for hyperelastic material modelling'
tags:
  - Julia
  - hyperelasticity
  - solid mechanics
authors:
  - name: Carson Farmer
    orcid: 0000-0000-0000-0000
    equal-contrib: false
    affiliation: 1
  - name: Hector Medina
    orcid: 0000-0000-0000-0000
    equal-contrib: false # (This is how you can denote equal contributions between multiple authors)
affiliations:
 - name: School of Engineering, Liberty University, Lynchburg, VA, United States
   index: 1
date: 28 March 2023
bibliography: paper.bib
---

# Summary

`Hyperelastics.jl` is a Julia [@Bezanson2017]  implementation for the largest (70+) collection of hyperelastic material models. The package provides a set of analytical and data-driven strain energy density functions (SEDF) and the tools required to calibrate the models to material tests. The package is designed to leverage multiple-dispatch to define a common set of functions for calculating the SEDF, Second Piola Kirchoff stress tensor, and the Cauchy stress tensor. The package provides: 1) a material model library that is AD compatible and 2) a set of extensible methods for easily defining and testing new material models. The package leverages the `ContinuumMechanicsBase.jl` pacakge for defining the continuum scale quantities and their corresponding relationships.

# Statement of need

The development of `Hyperelastics.jl` began as a study of the accuracy for a variety of material models for a set of experimental data. Often, researcher rely on custom implementations of material models and the data fitting process to find material parameters that match their experimental data. The SEDFs included in the package cover most of the available analytical models from the literature to date. Furthermore, a selection of data-driven models are incldued as a starting point for the development of new methods. 


`Hyperelastics.jl` is part of the Multi-Scale Material Modelling ($M^3$) Suite being developed in the Translational Robotics and Controls Engineering Research (TRACER) Lab at Liberty University. A pure Julia implementation allows for the use of automatic differentiation (AD) packages to calculate the partial derivatives of the SEDF. `Hyperelastics.jl` is designed to leverage multiple-dispatch to define a common set of functions for calculating the SED, Second Piola Kirchoff Stress Tensor, and the Cauchy Stress Tensor. The package provides a set of hyperelastic models and an interface to `Optimization.jl` for fitting model parameters.

# Availability

Hyperelastics.jl can be found on [github]

# Acknowledgements

The TRACER Lab is supported by the School of Engineering and the Center for Engineering Research and Education (CERE) at Liberty University.

# References
