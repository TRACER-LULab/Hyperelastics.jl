```@meta
CurrentModule = Hyperelastics
```

# API

```@index
```

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

# Helper Functions

```@autodocs
Modules = [Hyperelastics]
Order = [:function]
```

```@bibliography
```