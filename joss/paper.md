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

The development of `Hyperelastics.jl` began as a study of the accuracy for a variety of material models for a set of experimental data. Often, researcher rely on custom implementations of material models and the data fitting process to find material parameters that match their experimental data. 


`Hyperelastics.jl` is part of the Multi-Scale Material Modelling ($M^3$) Suite being developed in the Translational Robotics and Controls Engineering Research (TRACER) Lab at Liberty University. A pure Julia implementation allows for the use of automatic differentiation (AD) packages to calculate the partial derivatives of the SEDF. `Hyperelastics.jl` is designed to leverage multiple-dispatch to define a common set of functions for calculating the SED, Second Piola Kirchoff Stress Tensor, and the Cauchy Stress Tensor. The package provides: 1) a material model library that is AD compatible and 2) a set of extensible methods for easily defining and testing new material models.

# Functionality
The most basic definition in `Hyperelastics.jl` is the SEDF. The material models are implemented primarily by the SEDF with AD rules being defined to generate the Nominal, $S$, and True, $\sigma$ stresses felt in the material

$S_{ij} = 2\frac{\partial W}{\partial C_{ij}} = \frac{\partial W}{\partial \lambda_i}$

$\sigma_i = \frac{\lambda_i}{J}\frac{\partial W}{\partial \lambda_i}$
<!-- `Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`). -->

<!-- `Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike. -->



# Nomenclature:

| Variable | Term |
| :---: | :---: |
| $\psi$ | Strain energy density |
| $S$, $[\mathbf{S}]$ | Second Piola Kirchoff stress |
| $P$, $[\mathbf{P}]$ | First Piola Kirchoff stress |
| $\sigma$ | Cauchy stress |
| $\lambda_i$ | Principal stretch in direction $i$|
| $[\mathbf{F}]$ | Deformation gradient tensor |
| $[\mathbf{C}]$ | Right Cauchy-Green tensor |
| $[\mathbf{B}]$ | Left Cauchy-Green tensor |

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
<!-- - `@author:2001`  ->  "Author et al. (2001)" -->
<!-- - `[@author:2001]` -> "(Author et al., 2001)" -->
<!-- - `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)" -->

# Availability

Hyperelastics.jl can be found on [github]

# Figures

Figures can be included like this:
<!-- ![Caption for example figure.\label{fig:example}](figure.png) -->
<!-- and referenced from text using \autoref{fig:example}. -->

Figure sizes can be customized by adding an optional second parameter:
<!-- ![Caption for example figure.](figure.png){ width=20% } -->

# Acknowledgements

The TRACER Lab is supported by the School of Engineering and the Center for Engineering Research and Education (CERE) at Liberty University.

# References
