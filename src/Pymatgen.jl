module Pymatgen
    using PythonCall
    using ..Types


    export py_Struc, py_Incar
    export get_py_struc, get_sym_op_vec


    const py_Struc = PythonCall.pynew()
    const py_Incar = PythonCall.pynew()
    const py_Sga = PythonCall.pynew()
    
    function __init__()
        PythonCall.pycopy!(
            py_Struc,
            pyimport("pymatgen.core").Structure
        )
        PythonCall.pycopy!(
            py_Incar,
            pyimport("pymatgen.io.vasp.inputs").Incar
        )
        PythonCall.pycopy!(
            py_Sga,
            pyimport("pymatgen.symmetry.analyzer").SpacegroupAnalyzer
        )
    end


    get_py_struc(filepath) = py_Struc.from_file(filepath)


    function get_sym_op_vec(py_struc; symprec=1e-2, angle_tolerance=5.0)
        py_sga = py_Sga(
            py_struc,
            symprec=symprec,
            angle_tolerance=angle_tolerance
        )
        
        py_sym_dict = py_sga.get_symmetry_dataset()

        spg_num = pyconvert(Int64, py_sym_dict["number"])
        lattice_mat = permutedims(
            pyconvert(Array{Float64}, py_struc.lattice.matrix),
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