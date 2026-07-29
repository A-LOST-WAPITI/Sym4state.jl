# Sym4state

_Un programa diseñado específicamente para simplificar el cálculo de la matriz de interacción magnética y simular texturas de espín bajo diversas condiciones ambientales._

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://a-lost-wapiti.github.io/Sym4state.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://a-lost-wapiti.github.io/Sym4state.jl/dev/)
[![Build Status](https://github.com/a-lost-wapiti/Sym4state.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/a-lost-wapiti/Sym4state.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/a-lost-wapiti/Sym4state.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/a-lost-wapiti/Sym4state.jl)

## Instalación

Por favor, proceda a adquirir e instalar la versión apropiada de Julia (versión 1.9 o superior) para su plataforma específica visitando [esta ubicación](https://julialang.org/downloads/).

Y añada el paquete mediante:

```julia
julia> import Pkg

julia> Pkg.add(url="https://github.com/A-LOST-WAPITI/Sym4state.jl")
```

## Breve descripción general

Las funciones principales son `pre_process`, `post_process` y `Sym4state.MC.mcmc`.

La función `pre_process` facilita el cálculo optimizado de las matrices de interacción dentro de un radio de corte especificado. Para ilustrar su uso, consideremos un ejemplo donde queremos evaluar todas las interacciones dentro de un radio de corte de 5 Å para una monocapa de $\ce{CrI3}$. Para lograrlo, puede ejecutar sin esfuerzo el siguiente fragmento de código, asumiendo que tiene los archivos requeridos como POSCAR (que contiene la estructura de $\ce{CrI3}$), así como los archivos POTCAR, INCAR y KPOINTS utilizados en la determinación de la energía:

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

Tenga en cuenta que debe reemplazar `"POSCAR"`, `"POTCAR"`, `"INCAR"` y `"KPOINTS"` por los nombres de archivo reales correspondientes a los archivos de entrada específicos de su sistema.

Al invocar la función `pre_process`, esta automatizará la creación de configuraciones magnéticas con el número mínimo de energías requeridas para la determinación de la interacción. Estas configuraciones se guardarán en directorios separados. Una vez que todos los cálculos de DFT hayan convergido, puede emplear la función `post_process` para obtener la matriz de interacción. Este proceso puede realizarse utilizando el siguiente fragmento de código de Julia:

```julia
pair_mat, coeff_array = Sym4state.post_process("./cal.jld2")
```

La variable `pair_mat` contiene los índices de los pares de átomos, mientras que `coeff_array` contiene la matriz de interacción. Con esta información, puede establecer una simulación de Monte Carlo para analizar más a fondo el sistema mediante:

```julia
using Unitful, UnitfulAtomic, CUDA

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

Después de la simulación de Monte Carlo, se puede ilustrar la configuración magnética final simplemente mediante:

```julia
using CairoMakie
fig = Sym4state.MC.plot(
    states_over_env[:, :, :, :, end],
    Float32[1 -1/2; 0 sqrt(3)/2],
    Float32[1/3 2/3; 2/3 1/3],
    1:24,
    1:24
)
```

Es importante señalar que la matriz de celda 2D, las coordenadas fraccionales de los átomos magnéticos y el rango de trazado deseado de celdas a lo largo de los ejes x e y deben suministrarse como entradas a esta función.
