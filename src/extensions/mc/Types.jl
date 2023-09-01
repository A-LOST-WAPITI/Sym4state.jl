module MCTypes
    export Lattice, Environment, MCMethod

    mutable struct Lattice{T<:AbstractFloat}
        size::NTuple{2, Int}
        cell_mat::AbstractMatrix{T}
        offset_array::AbstractMatrix{T}
        magmom_vector::AbstractVector{T}
        point_idx_array::AbstractArray{Int, 3}
        interact_coeff_array::AbstractArray{T, 4}
    end

    struct Environment{T<:AbstractFloat}
        temperature::T
        magnetic_field::AbstractVector{T}
    end

    struct MCMethod
        step_equilibration_num::Int
        step_measure_num::Int
        #TODO: PTMC
    end
end