# Sym4state

_A program specifically designed to simplify the computation of magnetic interaction matrix and simulate spin textures under various environmental conditions._

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://a-lost-wapiti.github.io/Sym4state.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://a-lost-wapiti.github.io/Sym4state.jl/dev/)
[![Build Status](https://github.com/a-lost-wapiti/Sym4state.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/a-lost-wapiti/Sym4state.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/a-lost-wapiti/Sym4state.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/a-lost-wapiti/Sym4state.jl)

## Installation

Please proceed to acquire and install the appropriate version of Julia (version 1.9 or higher) for your specific platform by visiting [this location](https://julialang.org/downloads/).

And add the package by:

```julia
julia> import Pkg

julia> Pkg.add(url="https://github.com/A-LOST-WAPITI/Sym4state.jl")
```

## Brief overview

The main functions are `pre_process`, `post_process` and `Sym4state.MC.mcmc`.

The `pre_process` function facilitates the streamlined computation of interacting matrices within a specified cutoff radius. To illustrate its usage, let's consider an example where we aim to evaluate all the interactions within a cutoff radius of 5 Å for a monolayer of $\ce{CrI3}$. To achieve this, you can effortlessly execute the following code snippet, assuming you have the required files such as POSCAR (containing the structure of $\ce{CrI3}$), as well as the POTCAR, INCAR, and KPOINTS files utilized in energy determination:

```julia
using Sym4state
Sym4state.pre_process(
    "./POSCAR",
    [24],
    5.0;
    incar_path="./INCAR",
    potcar_path="./POTCAR",
    kpoints_path="./KPOINTS"
)
```

Note that you should replace `"POSCAR"`, `"POTCAR"`, `"INCAR"`, and `"KPOINTS"` with the actual filenames corresponding to your system's specific input files.

Upon invoking the `pre_process` function, it will automate the creation of magnetic configurations with the minimal number of energies required for interaction determination. These configurations will be saved in separate directories. Once all the DFT calculations have converged, you can employ the `post_process` function to obtain the interacting matrix. This process can be accomplished using the following Julia code snippet:

```julia
pair_mat, coeff_array = Sym4state.post_process("./cal.jld2")
```

The variable `pair_mat` contains the indices of atom pairs, while `coeff_array` holds the interacting matrix. With this information, you can establish a Monte Carlo simulation to further analyze the system by:

```julia
using Unitful, UnitfulAtomic, CUDA

mcconfig = Sym4state.MC.MCConfig{Float32}(
    lattice_size=[128, 128],
    magmom_vector=[3.5, 3.5],
    pair_mat=pair_mat,
    interact_coeff_array=coeff_array,
    temperature=collect(150:-2:0) * u"K" .|> austrip,
    magnetic_field=zeros(3),
    equilibration_step_num=100_000,
    measuring_step_num=100_000
)
(
    states_over_env,
    norm_mean_mag_over_env,
    susceptibility_over_env,
    specific_heat_over_env
) = Sym4state.MC.mcmc(
    mcconfig,
    backend=CUDABackend()
)
```