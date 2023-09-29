```@meta
CurrentModule = Hyperelastics
```

# API Reference

# Incompressible Models

```@autodocs
Modules = [Hyperelastics]
Filter = x -> typeof(x) === UnionAll && x <:Hyperelastics.AbstractIncompressibleModel
```

# Compressible Models

```@autodocs
Modules = [Hyperelastics]
Filter = x -> typeof(x) === UnionAll && x <:Hyperelastics.AbstractCompressibleModel
```

# Data Driven Models

```@autodocs
Modules = [Hyperelastics]
Filter = x -> typeof(x) === UnionAll && x <:Hyperelastics.AbstractDataDrivenHyperelasticModel
```

# Functions

```@autodocs
Modules = [Hyperelastics]
Order = [:function]
```