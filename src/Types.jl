"""
This module defines structures containing crystal structures,
symmetric operations as well as relation among different magnetic configurations.

This module exports:

$(EXPORTS)

"""
module Types

using DocStringExtensions
using Base
using StaticArrays
using LinearAlgebra
using PeriodicTable
using Printf
using CellListMap
using DataStructures: IntDisjointSets, find_root!


export Struc, SymOp, Map, CoeffMatRef

const A_IDX_MAT::Matrix{Vector{Int}} = [
    [[ 1,  2,  3,  4]] [[ 5,  6,  7,  8]] [[ 9, 10, 11, 12]];
    [[ 5,  6,  7,  8]] [[13,  2,  4, 14]] [[15, 16, 17, 18]];
    [[ 9, 10, 11, 12]] [[15, 16, 17, 18]] [[19,  1,  3, 20]]
]

sprintf(fmt::String, args...) = @eval Types @sprintf($(fmt), $(args...))

function vec2str(x; fmt::String="%7.4f")
    return "[" * join(
        [sprintf(fmt, item) for item in x],
        ", "
    ) * "]"
end

function mat2str(x; fmt::String="%10.4f")
    str = ""
    for row in eachrow(x)
        str *= join(
            [sprintf(fmt, item) for item in row],
            " "
        ) * "\n"
    end

    return str
end

"""
The data structure designed to store an atom includes its atomic number, position, and spin.

Two `Atom`s can be deemed approximately equal using the `isapprox` function if they share the same atomic number, possess position vectors that are approximately equal, and have spin vectors that are approximately equal.

$(TYPEDEF)

---

$(TYPEDFIELDS)
"""
struct Atom{T<:AbstractFloat}
    "Atomic number."
    num::Int
    "Position vector in fractional coordinates."
    pos::SVector{3, T}
    "Spin direction vector in cartesian coordinates."
    spin::SVector{3, T}
end

function Base.isapprox(x::Atom{T}, y::Atom{T}; check_spin::Bool=true, atol=T(1e-2)) where T
    num_flag = (x.num == y.num)
    spin_flag = !check_spin || isapprox(
        x.spin,
        y.spin,
        atol=atol
    )

    dis_1d_vec::MVector = @. abs(x.pos - y.pos)
    cir_flag_vec = dis_1d_vec .> 0.5
    dis_1d_vec[cir_flag_vec] .-= 1
    dis = norm(dis_1d_vec)
    pos_flag = dis < atol

    approx_flag = (num_flag && spin_flag && pos_flag)

    return approx_flag
end

function Base.show(io::IO, ::MIME"text/plain", atom::Atom)
    print(
        io,
        @sprintf("%3s", elements[atom.num].symbol) *
            " @ " * vec2str(atom.pos) *
            " with spin " * vec2str(atom.spin, fmt="%2d")
    )

    return nothing
end
Base.show(io::IO, atom::Atom) = show(io, "text/plain", atom)

"""
The data structure created for storing a material structure.

$(TYPEDEF)

---

$(TYPEDFIELDS)
"""
struct Struc
    "Unique index number to identify the material."
    uni_num::Int
    "The lattice matrix [a b c]."
    lattice_mat::SMatrix
    "Number of atoms in this structure."
    atom_count::Int
    "Atomic numbers for different atoms."
    num_vec::SVector
    "Fractional coordinates for different atoms."
    pos_mat::SMatrix
    "Spin directions for different atoms."
    spin_mat::SMatrix
end
"""
A convenient approach to construct a new [`Struc`](@ref) with essential parameters exclusively:

$(TYPEDSIGNATURES)
"""
function Struc(
    uni_num::Int,
    lattice_mat::AbstractMatrix{T},
    num_vec::AbstractVector{Int},
    pos_mat::AbstractMatrix{T},
    spin_mat::AbstractMatrix{T}
) where T
    atom_count = length(num_vec)
    @assert size(lattice_mat) == (3, 3)
    @assert size(pos_mat) == (3, atom_count)
    @assert size(spin_mat) == (3, atom_count)

    num_vec = SVector{atom_count, Int}(num_vec)
    lattice_mat = SMatrix{3, 3, T}(lattice_mat)
    pos_mat = SMatrix{3, atom_count, T}(pos_mat)
    spin_mat = SMatrix{3, atom_count, T}(spin_mat)

    return Struc(
        uni_num,
        lattice_mat,
        atom_count,
        num_vec,
        pos_mat,
        spin_mat
    )
end

function Base.iterate(struc::Struc, state=1)
    if state > struc.atom_count
        return nothing
    else
        return (
            struc[state],
            state + 1
        )
    end
end

Base.length(struc::Struc) = length(struc.atom_count)

Base.getindex(struc::Struc, idx::Int) = @views Atom(
    struc.num_vec[idx],
    struc.pos_mat[:, idx],
    struc.spin_mat[:, idx]
)
Base.getindex(struc::Struc, idics::Union{AbstractRange, AbstractVector}) = [
    struc[idx]
    for idx in idics
]

function Base.show(io::IO, ::MIME"text/plain", struc::Struc)
    println(io, "Structure Summary:")
    ## Lattice matrix
    print(
        io,
        "Lattice [a, b, c]\n" *
            mat2str(struc.lattice_mat)
    )
    ## Atoms
    println(
        io,
        "Atoms"
    )
    for atom in struc
        print(io, "  ")
        show(io, atom)
        println(io)
    end
end
Base.show(io::IO, struc::Struc) = show(io, "text/plain", struc)

"""
The data structure created for storing a symmetric operation.

$(TYPEDEF)

---

$(TYPEDFIELDS)
"""
struct SymOp{T<:AbstractFloat}
    "The rotational matrix for atoms."
    rot_mat::SMatrix{3, 3, T}
    "The rotational matrix for spins."
    spin_rot_mat::SMatrix{3, 3, T}
    "The translation vector for atoms."
    trans_vec::SVector{3, T}
    "The indicator flag to denote the propriety of the rotation."
    proper::Bool
    "The indicator flag to denote whether only translation operation is performed"
    trans_only::Bool
    "The indicator flag to denote the presence of a time reversal operation."
    time_rev::Bool
end
"""
A convenient approach to construct a new [`SymOp`](@ref) with essential parameters exclusively:

$(TYPEDSIGNATURES)

---

As the spins are pseudovectors, their rotational matrix can be defined as:
```math
\\mathbf{R}^\\prime = T \\cdot \\|\\mathbf{R}\\| \\cdot \\mathbf{R}
```
where ``\\mathbf{R}`` is the rotational matrix for atoms and ``T`` is the time reversal operation.
"""
function SymOp(
    rot_mat::AbstractMatrix{T},
    trans_vec::AbstractVector{T},
    time_rev::Bool,
    lattice_mat::AbstractMatrix{T},
    inv_lattice_mat::AbstractMatrix{T}
) where T
    @assert size(rot_mat) == (3, 3) "The rotation matrix should have a shape of 3x3!"
    @assert length(trans_vec) == 3

    rot_mat = SMatrix{3, 3, T}(rot_mat)
    trans_vec = SVector{3, T}(trans_vec)
    proper::Bool = det(rot_mat) > 0
    spin_rot_mat = SMatrix{3, 3, T}(
        (
            time_rev ? -one(T) : one(T)
        ) * (
            proper ? one(T) : -one(T)
        ) * lattice_mat * rot_mat * inv_lattice_mat
    )

    return SymOp(
        rot_mat,
        spin_rot_mat,
        trans_vec,
        proper,
        trans_only,
        time_rev
    )
end

Base.:*(op::SymOp, struc::Struc) = begin
    new_pos_mat = mod.(op.rot_mat * struc.pos_mat .+ op.trans_vec, 1)
    new_spin = op.spin_rot_mat * struc.spin_mat

    return Struc(
        struc.uni_num,
        struc.lattice_mat,
        struc.atom_count,
        struc.num_vec,
        new_pos_mat,
        new_spin
    )
end

Base.:*(op::SymOp, mat::AbstractMatrix) = begin
    @assert size(mat) == (3, 3)
    new_mat = op.spin_rot_mat * mat * transpose(op.spin_rot_mat)

    return new_mat
end

function Base.show(io::IO, ::MIME"text/plain", sym_op::SymOp)
    function vec2str(x)
        return "[" * join(
            [@sprintf("%6.4f", item) for item in x],
            ", "
        ) * "]"
    end
    function mat2str(x)
        str = ""
        for row in eachrow(x)
            str *= join(
                [@sprintf("%10.4f", item) for item in row],
                " "
            ) * "\n"
        end

        return str
    end

    println(io, "Symmtry operation summary:")
    ## Rotation matrix
    proper = sym_op.proper > 0 ? "Proper" : "Improper"
    println(
        io,
        "Rotation matrix ($(proper)):\n" *
            mat2str(sym_op.rot_mat)
    )
    ## Translation vector
    println(
        io,
        "Translation vector: " *
            vec2str(sym_op.trans_vec)
    )
    ## Time reversal
    time_rev = sym_op.time_rev ? "True" : "False"
    println(
        io,
        "Time reversal: " * time_rev
    )

    return nothing
end
Base.show(io::IO, sym_op::SymOp) = show(io, "text/plain", sym_op)

"""
The data structure created for storing all the configurations needed to calculate an interaction matrix.

$(TYPEDEF)

---

$(TYPEDFIELDS)
"""
struct Map
    "A matrix encompassing all the `uni_num` values of the configurations required to compute the elements of the interaction matrix."
    map_mat::Matrix{Vector{Int8}}
    "A vector that stores the connections between the indices of the required real configurations within the set of all configurations."
    fallback_vec::Vector{Int8}
    "A vector that stores all the magnetic configurations."
    struc_vec::Vector{Struc}
    "The type of this interaction matrix. 1 for exchange interaction matrix while 2 for single-ion anisotropy matrix."
    type::Int
end
"""
A convenient approach to construct a new [`Map`](@ref) with essential parameters exclusively:

$(TYPEDSIGNATURES)
"""
function Map(
    fallback_ds::IntDisjointSets,
    all_struc_vec::Vector{Struc};
    rotation_symmetry_flag=false
)
    parents = fallback_ds.parents
    if length(parents) == 36
        type = 1    # J matrix
        temp = eachslice(
            reshape(
                parents,
                (4, 3, 3)
            ),
            dims=(2, 3)
        )
    else length(parents) == 20
        type = 2    # A matrix
        temp = deepcopy(A_IDX_MAT)
        for (element_idx, element_comp) in enumerate(temp)
            if rotation_symmetry_flag && element_idx != 9 # only Azz - Axx is nonzero
                coeff = 0
            else
                coeff = 1
            end

            @. element_comp = coeff * parents[element_comp]
        end
    end

    map_mat = [zeros(Int8, 4) for _ = 1:3, _ = 1:3]
    for idx in eachindex(temp)
        element_comp = temp[idx]
        part1 = element_comp[[1, 4]]
        part2 = element_comp[[2, 3]]

        part_diff = setdiff(part1, part2)
        if length(part_diff) != 0
            map_mat[idx] .= element_comp
        end
    end

    fallback_vec = filter(>(0), vcat(map_mat...)) |> unique |> sort

    return Map(
        map_mat,
        fallback_vec,
        all_struc_vec[fallback_vec],
        type
    )
end

function Base.show(io::IO, ::MIME"text/plain", map::Map)
    for row in eachrow(map.map_mat)
        for item in row
            print(io, vec2str(item, fmt="%2d"))
            print(io, "\t")
        end
        println(io)
    end
    type_str = isone(map.type) ? "J" : "A"

    println(io, "A reduced \"$(type_str)\" matrix with $(length(map.fallback_vec)) unique configurations.")
end
Base.show(io::IO, map::Map) = show(io, "text/plain", map)

"""
The data structure created for storing the equivalent relation between two atom pairs.

$(TYPEDEF)

---

$(TYPEDFIELDS)
"""
struct CoeffMatRef
    "Unique index number to identify the group of those two atom pairs."
    group_idx::Int
    "A vector indicating two atom pairs."
    pair_vec::AbstractVector{Int}
    "The operation linking those two atom pairs."
    op::SymOp
end


end