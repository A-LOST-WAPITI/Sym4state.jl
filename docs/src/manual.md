# Manual

```@meta
CurrentModule = Sym4state.ModCore
```

```@setup pre_and_post
using PrintFileTree
local pair_mat, coeff_array
```

For the theoretical exploration of the magnetic properties of magnets, the bilinear Heisenberg model proves to be a useful framework for representing magnetic interactions, which can be described by

```math
\mathcal{H} = \sum_{i < j} S_i \cdot \mathcal{J}_{i j} \cdot S_j + \sum_{i} S_i \cdot \mathcal{A} \cdot S_i - m \sum_{i} S_i \cdot \vec{B}
```

where the symbol ``\mathcal{J}_{ij}`` denotes the exchange interaction matrix between two spins, ``S_i`` and ``S_j``, the matrix ``\mathcal{A}`` represents the single-ion anisotropy. To determine the magnetic interaction matrix elements, researchers often employ the four-state method [^1] [^2] [^3]. This method involves calculating the energies of four distinct magnetic configurations, allowing the extraction of individual components for the exchange matrix.

Extending this method to each element of the exchange matrix requires calculating a total of 36 energies to obtain the complete matrix. It should be noted that some energies are degenerate due to the symmetry of the material. Nonetheless, performing a manual symmetric analysis to streamline the number of energy calculations remains a challenging endeavor, as there exists a potential risk of omitting or misinterpreting certain symmetric operations.

## Pre-process

One can use our program to streamline the simpilifing and calculating process easily. For example, with a POSCAR file of monolayer ``\ce{CrI3}``

```raw
Cr2 I6                                  
   1.00000000000000     
     7.1131374882967124    0.0000000000000000    0.0000000000000000
    -3.5565687441483571    6.1601577654763897    0.0000000000000000
     0.0000000000000000    0.0000000000000000   18.0635365764484419
   Cr   I 
     2     6
Direct
  0.6666666666666643  0.3333333333333357  0.5000000247180765
  0.3333333333333357  0.6666666666666643  0.5000000501683317
  0.6415738047516142  0.9999977877949036  0.4116659127023310
  0.3584239830432894  0.3584261952483858  0.4116659127023310
  0.0000022122051035  0.6415760169567106  0.4116659127023310
  0.3584241488090230  0.9999980859273947  0.5883340783387269
  0.6415739371183646  0.6415758511909699  0.5883340783387269
  0.0000019140726053  0.3584260628816354  0.5883340783387269
```

and the proper setted INCAR, POTCAR and KPOINTS for making SCF calculation, one can simply using `Sym4state.jl` to generate all the input files for calculating the nearest exchange interaction and the single-ion anisotropy interaction as follows:

```@example pre_and_post
using Sym4state
cd("CrI3") do   # hide
Sym4state.pre_process(
    "./POSCAR",
    [24],   # Take Cr element as magnetic
    5.0     # There exists an interaction between atoms within a distance of 5 Å.
)
end # hide
```

This function will utilize the [`supercell_check`](@ref) method to create a supercell for the provided structure. The supercell should be sufficiently large to ensure that no more than one connection exists within a specified cutoff radius between any two atoms. For the given case of a monolayer of ``\ce{CrI3}`` with a cutoff radius of 5 Å, a ``2 \times 2 \times 1`` supercell will provide sufficient size. The supercell diagram below labels all the ``\ce{Cr}`` atoms:

![Top view of monolayer ``\ce{CrI3}``](figs/CONTCAR.webp)

Within the 5 Å cutoff radius, the monolayer of ``\ce{CrI3}`` exhibits two distinct groups of interactions. The first group corresponds to interactions between nearest neighbors, whereas the second group pertains to interactions arising from single-ion anisotropy. It is important to note that all atom pairs within the same group are considered equivalent. This equivalence implies the existence of symmetric operations that can transform one interaction matrix into another, highlighting the underlying symmetry of the system.

As evidenced by the output obtained from the `pre_process` function, the initial group contains 6 pairs that are equivalent, while the second group consists of 2 equivalent pairs. Despite the potential for simplifying the calculations involving various interaction matrices through the use of symmetric operations, there remains one particular interaction matrix that necessitates the calculation of the fewest number of configurations. In the case of the nearest neighbor interaction, it is essential to compute the energies for a minimum of 9 magnetic configurations. Conversely, when dealing with the single-ion anisotropy interaction, the energies of at least 2 magnetic configurations need to be evaluated.

The function will restore all the relations between different energies and configurations into a file `cal.jld2`. Moreover, this function will generate numerous directories to store input files corresponding to the various magnetic configurations.

```@example pre_and_post
printfiletree("CrI3")   # hide
```

All the path of those directories is stored in the file `cal_list`, one could use this file to create a [Slurm](https://slurm.schedmd.com/)'s job array by submitting a shell like:

```bash
#!/bin/sh

#SBATCH -n 144
#SBATCH --array=1-11%2

module load vasp-6.3.2-optcell

target_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" cal_dir_list)

cd ${target_dir}

srun vasp_ncl
```

This shell script aims to create a Slurm job array to compute the energies of all 11 magnetic configurations, while efficiently managing computational resources by allowing a maximum of 2 jobs to run simultaneously.

## Post-process

Once all the calculations have converged, you can utilize the `post_process` function to extract the energies associated with different configurations. This process ultimately leads to the construction of an interaction matrix.

```@example pre_and_post
cd("CrI3") do   # hide
global pair_mat, coeff_array    # hide
mv("../oszicar.tar.gz", "./oszicar.tar.gz") # hide
run(`tar -xvzf oszicar.tar.gz`) # hide
for (idx, dir_name) in enumerate(readlines("cal_dir_list")) # hide
    cp("oszicar/OSZICAR_$(idx)", dir_name * "OSZICAR")  # hide
end # hide
pair_mat, coeff_array = Sym4state.post_process("./cal.jld2")
end # hide
```

We can examine the dimensions of `pair_mat` and `coeff_array`, which store the indices of the starting and ending points for various atom pairs and their corresponding interaction matrices, respectively.

```@repl pre_and_post
size(pair_mat)
size(coeff_array)
```

Hence, we observe that there exist a total of 8 interactions within a cutoff radius of 5 Å. Let us inspect a specific entry in `pair_mat` that contains the indices representing an atom pair:

```@repl pre_and_post
pair_mat[:, 1]
```

The initial and final numbers correspond to the indices of the starting and ending point atoms, respectively. The second and third numbers indicate the offset of the primitive cell along the x-axis and y-axis.

## Monte Carlo Simulation

With the former result `pair_mat` and `coeff_array`, we could set up a configuration for Monte Carlo simulation to determining the phase transition temperature or magnetic texture like:

```@repl pre_and_post
using Unitful, UnitfulAtomic
mcconfig = Sym4state.MC.MCConfig{Float32}(
    lattice_size=[128, 128],
    magmom_vector=[3.5, 3.5],
    pair_mat=pair_mat,
    interact_coeff_array=coeff_array,
    temperature=collect(150:-2:0),
    magnetic_field=zeros(3),
    equilibration_step_num=100_000,
    measuring_step_num=100_000
)
```

In the aforementioned code snippet, we have configured a simulated annealing simulation, commencing at a temperature of 150 K and progressively reducing it to 0 K in steps of 2 K. The simulation operates on a ``128 \times 128`` supercell of ``\ce{CrI3}`` using the previously computed interaction matrix. To assess the system, we perform a preliminary equilibration phase consisting of ``100000`` sweeps, followed by a measurement phase comprising ``100000`` sweeps for acquiring physical quantities. It is worth noting that the magnetic field is absent, rendering the `magmom_vector` inconsequential.

With the created `mcconfig`, one can initiate a Monte Carlo simulation as follows:

```julia
(
    states_over_env,
    norm_mean_mag_over_env,
    susceptibility_over_env,
    specific_heat_over_env
) = Sym4state.MC.mcmc(
    mcconfig,
    backend=Sym4state.MC.CPU()
    progress_enabled=false,
    log_enabled=false
)
```

The parameter `backend` can be configured to employ `CUDABackend()` provided by [`CUDA.jl`](https://github.com/JuliaGPU/CUDA.jl) or any other backends supported by [`KernelAbstractions.jl`](https://github.com/JuliaGPU/KernelAbstractions.jl) to enhance performance utilizing the GPU.

The `MCConfig` can also be stored into a `.toml` file by:

```@example pre_and_post
cd("CrI3") do   # hide
Sym4state.MC.save_config("CrI3.toml", mcconfig)
end # hide
```

or it can also be restored by:

```@example pre_and_post
cd("CrI3") do   # hide
mcconfig = Sym4state.MC.load_config("CrI3.toml")
end # hide
```

## Functions

```@docs
reduce_interact_mat_for_a_pair
supercell_check
pre_process
post_process
```

```@eval
rm("CrI3", recursive=true)
nothing
```

[^1]: Xiang, H. J., et al. "Predicting the spin-lattice order of frustrated systems from first principles." Physical Review B 84.22 (2011): 224429.
[^2]: Šabani, D., C. Bacaksiz, and M. V. Milošević. "Ab initio methodology for magnetic exchange parameters: Generic four-state energy mapping onto a Heisenberg spin Hamiltonian." Physical Review B 102.1 (2020): 014457.
[^3]: Xiang, Hongjun, et al. "Magnetic properties and energy-mapping analysis." Dalton Transactions 42.4 (2013): 823-853.
