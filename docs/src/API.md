```@meta
CurrentModule = Hyperelastics
```

# API

```@index
```

# Incompressible Models

```@autodocs
Modules = [Hyperelastics]
Filter = x -> typeof(x) === DataType && x <: Hyperelastics.AbstractIncompressibleModel
```

# Compressible Models

```@autodocs
Modules = [Hyperelastics]
Filter = x -> typeof(x) === DataType && x <: Hyperelastics.AbstractCompressibleModel

```

# Helper Functions

```@autodocs
Modules = [Hyperelastics]
Order = [:function]
```

```@bibliography
```