# Utilities

```@meta
CurrentModule = Sym4state.Utils
```

This page shows all the functions defined in `Sym4state.Utils`.

## Module

```@docs
Utils
```

## Functions

### Structures and configurations

This section contains functions working with structures or magnetic configurations.

```@docs
mag_config
py_struc_to_struc
get_all_interact_struc_vec
struc_compare
magonly
```

### Symmetric operations

This section contains functions working with space group and symmetric operations.

```@docs
check_z_rot_mat
check_z_rot
get_sym_op_vec
```

### Atom pairs

This section contains functions working with atom pairs.

```@docs
consider_pair_vec_in_radius
linear_idx_to_vec
get_corresponding_pair_vec
equal_pair
get_fixed_pair_vec
```

### Inputs and outputs

```@docs
to_vasp_inputs
get_py_struc
grep_energy
get_coeff_mat
get_one_interact_coeff_mat
get_all_interact_coeff_under_sym
get_pair_and_coeff
```

### Miscellaneous

```@docs
set_rwigs
check_unit_cell
```