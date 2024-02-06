# Manual

```@meta
CurrentModule = Sym4state.ModCore
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

```@example
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

![Top view of monolayer ``\ce{CrI3}``](figs/CONTCAR.png)

Within the 5 Å cutoff radius, the monolayer of ``\ce{CrI3}`` exhibits two distinct groups of interactions. The first group corresponds to interactions between nearest neighbors, whereas the second group pertains to interactions arising from single-ion anisotropy. It is important to note that all atom pairs within the same group are considered equivalent. This equivalence implies the existence of symmetric operations that can transform one interaction matrix into another, highlighting the underlying symmetry of the system.

The function above will create lots of directories storing input files for different magnetic configurations.

```@example
using PrintFileTree # hide
printfiletree("CrI3")   # hide
```

All the path of those directories is stored in the file `cal_list`, one could use this file to create a [Slurm](https://slurm.schedmd.com/)'s job array by submitting a shell like:

```bash
#!/bin/sh

#SBATCH -n 144
#SBATCH --array=1-2%2

module load vasp-6.3.2-optcell

target_dir=$(sed -n "${SLURM_ARRAY_TASK_ID}p" cal_dir_list)

cd ${target_dir}

srun vasp_ncl
```

## Post-process

## Module

```@docs
```

## Functions

```@docs
reduce_interact_mat_for_a_pair
supercell_check
```

```@eval
rm("CrI3", recursive=true)
nothing
```

[^1]: Xiang, H. J., et al. "Predicting the spin-lattice order of frustrated systems from first principles." Physical Review B 84.22 (2011): 224429.
[^2]: Šabani, D., C. Bacaksiz, and M. V. Milošević. "Ab initio methodology for magnetic exchange parameters: Generic four-state energy mapping onto a Heisenberg spin Hamiltonian." Physical Review B 102.1 (2020): 014457.
[^3]: Xiang, Hongjun, et al. "Magnetic properties and energy-mapping analysis." Dalton Transactions 42.4 (2013): 823-853.
