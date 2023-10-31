module MCTypes
    export Lattice, Environment, MCMethod

    mutable struct Lattice{T<:AbstractFloat}
        size::AbstractVector{Int}
        cell_mat::AbstractMatrix{T}
        offset_mat::AbstractMatrix{T}
        magmom_vector::AbstractVector{T}
        pair_mat::AbstractMatrix{Int}
        interact_coeff_array::AbstractArray{T, 3}
    end

    struct Environment{T<:AbstractFloat}
        temperature::T
        magnetic_field::AbstractVector{T}
    end

    struct MCMethod
        equilibration_step_num::Int
        measuring_step_num::Int
        #TODO: PTMC
    end
end