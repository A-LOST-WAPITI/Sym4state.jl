# Types

```@meta
CurrentModule = Sym4state.Types
```

This page shows all the types defined in `Sym4state.Types`.

## Module

```@docs
Types
```

## Types and methods

```@docs
Atom{T<:AbstractFloat}
```

```@docs
Struc
Struc(::Int)
```

```@docs
SymOp{T<:AbstractFloat}
SymOp(::AbstractMatrix)
```

```@docs
Map
Map(
    fallback_ds::IntDisjointSets,
    all_struc_vec::Vector{Struc};
    rotation_symmetry_flag=false
)
```

```@docs
CoeffMatRef
```