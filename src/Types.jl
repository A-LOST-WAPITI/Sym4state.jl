module Types
    using Base
    using StaticArrays
    using LinearAlgebra
    using PeriodicTable
    using Printf
    using CellListMap
    using DataStructures: IntDisjointSets, find_root!


    export Struc, SymOp, FallbackList, Map, CoeffMatRef

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


    struct Atom
        num::Int64
        pos::SVector{3, Float64}
        spin::SVector{3, Float64}
    end

    Base.isapprox(x::Atom, y::Atom; atol=1e-2) = begin
        num_flag = (x.num == y.num)
        spin_flag = isapprox(
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

    struct Struc
        uni_num::Int64
        lattice_mat::SMatrix
        atom_count::Int64
        num_vec::SVector
        pos_mat::SMatrix
        spin_mat::SMatrix
    end
    function Struc(
        uni_num::Int64,
        lattice_mat::AbstractMatrix,
        num_vec::AbstractVector,
        pos_mat::AbstractMatrix,
        spin_mat::AbstractMatrix
    )
        atom_count = length(num_vec)
        @assert size(lattice_mat) == (3, 3)
        @assert size(pos_mat) == (3, atom_count)
        @assert size(spin_mat) == (3, atom_count)

        num_vec = SVector{atom_count, Int64}(num_vec)
        lattice_mat = SMatrix{3, 3, Float64}(lattice_mat)
        pos_mat = SMatrix{3, atom_count, Float64}(pos_mat)
        spin_mat = SMatrix{3, atom_count, Float64}(spin_mat)

        return Struc(
            uni_num,
            lattice_mat,
            atom_count,
            num_vec,
            pos_mat,
            spin_mat
        )
    end

    Base.iterate(struc::Struc, state=1) = begin
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
        @views Atom(
            struc.num_vec[idx],
            struc.pos_mat[:, idx],
            struc.spin_mat[:, idx]
        )
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


    struct SymOp
        rot_mat::SMatrix{3, 3, Float64}
        spin_rot_mat::SMatrix{3, 3, Float64}
        trans_vec::SVector{3, Float64}
        proper::Int8
        time_rev::Int8
    end
    function SymOp(rot_mat, trans_vec, time_rev, lattice_mat::Matrix{Float64}, inv_lattice_mat::Matrix{Float64})
        @assert size(rot_mat) == (3, 3) "The rotation matrix should have a shape of 3x3!"
        @assert length(trans_vec) == 3

        rot_mat = SMatrix{3, 3, Float64}(rot_mat)
        trans_vec = SVector{3, Float64}(trans_vec)
        proper::Int8 = sign(det(rot_mat))
        time_rev::Int8 = time_rev
        spin_rot_mat = SMatrix{3, 3, Float64}(
            time_rev * proper * lattice_mat * rot_mat * inv_lattice_mat
        )

        return SymOp(
            rot_mat,
            spin_rot_mat,
            trans_vec,
            proper,
            time_rev
        )
    end

    Base.:*(op::SymOp, struc::Struc) = begin
        new_pos_mat = mod1.(op.rot_mat * struc.pos_mat .+ op.trans_vec, 1)
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
        time_rev = sym_op.time_rev > 0 ? "False" : "True"
        println(
            io,
            "Time reversal: " * time_rev
        )

        return nothing
    end
    Base.show(io::IO, sym_op::SymOp) = show(io, "text/plain", sym_op)


    struct Map
        map_mat::Matrix{Vector{Int8}}
        fallback_vec::Vector{Int8}
        struc_vec::Vector{Struc}
        type::Int
    end

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


    struct CoeffMatRef
        group_idx::Int
        pair_vec::AbstractVector{Int}
        op::SymOp
    end
end