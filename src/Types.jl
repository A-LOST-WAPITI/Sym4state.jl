module Types
    using Base
    using StaticArrays
    using LinearAlgebra
    using PeriodicTable
    using Printf


    export Struc, SymOp, FallbackList


    struct Atom
        num::Int64
        pos::SVector{3, Float64}
        spin::SVector{3, Float64}
    end

    Base.isapprox(x::Atom, y::Atom; atol=1e-2) = begin
        num_flag = (x.num == y.num)
        pos_flag = isapprox(
            x.pos,
            y.pos,
            atol=atol
        )
        spin_flag = isapprox(
            x.spin,
            y.spin,
            atol=atol
        )

        approx_flag = (num_flag && pos_flag && spin_flag)

        return approx_flag
    end

    function Base.show(io::IO, ::MIME"text/plain", atom::Atom)
        function vec2str(x)
            return "[" * join(
                [@sprintf("%7.4f", item) for item in x],
                " "
            ) * "]"
        end

        println(
            io,
            "$(elements[atom.num].symbol)" *
                " @ " * vec2str(atom.pos) *
                " with spin " * vec2str(atom.spin)
        )

        return nothing
    end
    Base.show(io::IO, atom::Atom) = show(io, "text/plain", atom)

    struct Struc
        uni_num::Int64
        lattice_mat::SMatrix
        inv_lattice_mat::SMatrix
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

        inv_lattice_mat = inv(lattice_mat)

        num_vec = SVector{atom_count, Int64}(num_vec)
        lattice_mat = SMatrix{3, 3, Float64}(lattice_mat)
        inv_lattice_mat = SMatrix{3, 3, Float64}(inv_lattice_mat)
        pos_mat = SMatrix{3, atom_count, Float64}(pos_mat)
        spin_mat = SMatrix{3, atom_count, Float64}(spin_mat)

        return Struc(
            uni_num,
            lattice_mat,
            inv_lattice_mat,
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
                Atom(
                    struc.num_vec[state],
                    struc.pos_mat[:, state],
                    struc.spin_mat[:, state]
                ),
                state + 1
            )
        end
    end

    Base.length(struc::Struc) = struc.atom_count

    Base.isapprox(x::Struc, y::Struc; atol=1e-2) = begin
        approx_flag = true

        for one_atom in x
            temp_vec = [
                isapprox(one_atom, another_atom, atol=atol)
                for another_atom in y
            ]

            if any(temp_vec)
            else
                approx_flag = false
                break
            end
        end

        return approx_flag
    end

    function Base.show(io::IO, ::MIME"text/plain", struc::Struc)
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
        end
    end
    Base.show(io::IO, struc::Struc) = show(io, "text/plain", struc)


    struct SymOp
        rot_mat::SMatrix{3, 3, Float64}
        trans_vec::SVector{3, Float64}
        proper::Int8
        time_rev::Int8
    end
    function SymOp(rot_mat, trans_vec, time_rev)
        @assert size(rot_mat) == (3, 3) "The rotation matrix should have a shape of 3x3!"
        @assert length(trans_vec) == 3

        rot_mat = SMatrix{3, 3, Float64}(rot_mat)
        trans_vec = SVector{3, Float64}(trans_vec)
        proper::Int8 = (det(rot_mat) > 0) ? 1 : -1
        time_rev::Int8 = time_rev

        return SymOp(
            rot_mat,
            trans_vec,
            proper,
            time_rev
        )
    end

    Base.:*(op::SymOp, struc::Struc) = begin
        new_pos_mat = op.rot_mat * struc.pos_mat .+ op.trans_vec
        frac_spin = struc.inv_lattice_mat * struc.spin_mat
        new_frac_spin = op.time_rev * op.proper * op.rot_mat * frac_spin
        new_spin = struc.lattice_mat * new_frac_spin

        return Struc(
            struc.uni_num,
            struc.lattice_mat,
            struc.inv_lattice_mat,
            struc.atom_count,
            struc.num_vec,
            new_pos_mat,
            new_spin
        )
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


    struct FallbackList
        len::Int64
        tree::MVector
    end
    FallbackList(len::Int64) = FallbackList(len, @MVector zeros(Int64, len))
    function (fallback::FallbackList)(idx::Int64)
        @assert idx <= fallback.len

        parent_idx = fallback.tree[idx]
        if iszero(parent_idx)
            return idx
        else
            fallback(parent_idx)
        end
    end

    function update!(fallback::FallbackList, idx, target_idx)
        @assert 0 < idx <= fallback.len
        @assert 0 < target_idx <= fallback.len

        target_idx_now = fallback.tree[idx]
        if iszero(target_idx_now)
            fallback.tree[idx] = target_idx
        else
            target_idx = min(target_idx, target_idx_now)
            fallback.tree[idx] = target_idx
        end
    end
end