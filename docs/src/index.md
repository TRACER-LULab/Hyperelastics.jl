```@meta
CurrentModule = Hyperelastics
```

# ContinuumMechanicsBase

`ContinuumMechanicsBase` contains the core functionality for continuum mechanics methods. The methods are intended to be overloaded by the package implementing the common interface. The objective of the package is to provide a common interface between interacting continuum modelling packages. The package is currently used by TRACER Lab for multi-scale modelling of soft materials. Packages which rely on ContinuumMechanicsBase include:

- Hyperelastics.jl

Further packages are in development and will be made public upon completion. Please submit a PR if you are using this package in your project. 

```@index
```

```@autodocs
Modules = [Hyperelastics]
```

```@bibliography
```