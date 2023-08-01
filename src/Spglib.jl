module Spglib
    using PythonCall
    using ..Types


    export get_sym_op_vec


    const py_spglib = PythonCall.pynew()
    
    function __init__()
        PythonCall.pycopy!(
            py_spglib,
            pyimport("spglib")
        )
    end


    function get_sym_op_vec(py_struc; symprec=1e-2)
        py_lattice_mat = py_struc.lattice.matrix
        py_pos_mat = py_struc.frac_coords
        py_num_vec = py_struc.atomic_numbers
        cell = pytuple((py_lattice_mat, py_pos_mat, py_num_vec))
        
        py_sym_dict = py_spglib.get_symmetry_dataset(cell, symprec=symprec)

        spg_num = pyconvert(Int64, py_sym_dict["number"])
        lattice_mat = permutedims(
            pyconvert(Array{Float64}, py_lattice_mat),
            (2, 1)
        )
        inv_lattive_mat = inv(lattice_mat)
        rot_array = permutedims(
            pyconvert(Array{Float64}, py_sym_dict["rotations"]),
            (2, 3, 1)
        )
        trans_mat = permutedims(
            pyconvert(Array{Float64}, py_sym_dict["translations"]),
            (2, 1)
        )

        op_vec = SymOp[]
        num_op = size(rot_array, 3)
        for idx = 1:num_op
            rot_mat = rot_array[:, :, idx]
            trans_vec = trans_mat[:, idx]

            for time_rev = 1:-2:-1
                op = SymOp(
                    rot_mat,
                    trans_vec,
                    time_rev,
                    lattice_mat,
                    inv_lattive_mat
                )

                push!(op_vec, op)
            end
        end

        return spg_num, op_vec
    end
end