module MCTypes


using Parameters
using ArgCheck
using Unitful: @u_str, ustrip
using UnitfulAtomic: auconvert, austrip


export MCConfig


const MU_B::Rational{Int} = 1//2


@with_kw struct MCConfig{T<:Real}
    # parameters for lattice
    lattice_size::Vector{Int}   = [128, 128]
    magmom_vector::Vector{T}
    pair_mat::Matrix{Int}
    interact_coeff_array::Array{T, 3}
    # parameters for environment
    temperature::Vector{T}      = zeros(eltype(magmom_vector), 1)
    magnetic_field::Array{T} = zeros(eltype(magmom_vector), 3)
    # parameters for monte carlo
    equilibration_step_num::Int = 100_000
    measuring_step_num::Int     = 10_000
    decorrelation_step_num::Int = 10

    # dimension check when construct
    MCConfig{T}(
        lattice_size,
        magmom_vector,
        pair_mat,
        interact_coeff_array,
        temperature,
        magnetic_field,
        equilibration_step_num,
        measuring_step_num,
        decorrelation_step_num
    ) where {T} = begin
        @argcheck size(pair_mat, 2) == size(interact_coeff_array, 3) DimensionMismatch

        # turn parameters into atomic unit
        au_interact_coeff_array = @. austrip(auconvert(interact_coeff_array * u"meV"))
        au_temperature = @. austrip(auconvert(temperature * u"K"))
        au_magnetic_field = @. austrip(auconvert(MU_B * magnetic_field * u"T"))

        new(
            lattice_size,
            magmom_vector,
            pair_mat,
            au_interact_coeff_array,
            au_temperature,
            au_magnetic_field,
            equilibration_step_num,
            measuring_step_num,
            decorrelation_step_num
        )
    end

    MCConfig{T}(
        magmom_vector,
        pair_mat,
        interact_coeff_array
    ) where {T} = MCConfig{T}(
        magmom_vector=magmom_vector,
        pair_mat=pair_mat,
        interact_coeff_array=interact_coeff_array
    )
end
# Outer constructor for not given type specifically
MCConfig(
    magmom_vector::Vector{T},
    pair_mat::Matrix{Int},
    interact_coeff_array::Array{T, 3}
) where {T} = MCConfig{T}(
    magmom_vector=magmom_vector,
    pair_mat=pair_mat,
    interact_coeff_array=interact_coeff_array
)


end