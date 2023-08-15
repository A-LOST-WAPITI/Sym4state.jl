module Utils
    using InvertedIndices
    using LinearAlgebra
    using DelimitedFiles
    using PythonCall
    using ..Types
    using ..Pymatgen

    
    export fourstate, to_vasp_inputs


    const py_np = PythonCall.pynew()

    function __init__()
        PythonCall.pycopy!(
            py_np,
            pyimport("numpy")
        )
    end


    function mag_config(mag_count, target_idx_vec)
        @assert length(target_idx_vec) == 2
        @assert target_idx_vec[1] != target_idx_vec[2]

        axes_vec = [1, 2, 3]
        mag_config_array = zeros(3, mag_count, 36)
        config_count = 0
        for alpha in axes_vec, beta in axes_vec
            left_axes = setdiff(axes_vec, alpha, beta)[end]

            for sign_1 = 1:-2:-1, sign_2 = 1:-2:-1
                config_count += 1
                mag_alpha = zeros(3)
                mag_beta = zeros(3)
                mag_left = zeros(3)

                mag_alpha[alpha] = sign_1
                mag_beta[beta] = sign_2
                mag_left[left_axes] = 1

                mag_config_array[:, target_idx_vec[1], config_count] .= mag_alpha
                mag_config_array[:, target_idx_vec[2], config_count] .= mag_beta
                mag_config_array[:, Not(target_idx_vec), config_count] .= mag_left
            end
        end

        return mag_config_array
    end


    function fourstate(py_struc, mag_num_vec, target_idx_vec)
        py_lattice_mat = py_struc.lattice.matrix
        py_pos_mat = py_struc.frac_coords
        py_num_vec = py_struc.atomic_numbers

        lattice_mat = permutedims(
            pyconvert(Matrix{Float64}, py_lattice_mat),
            (2, 1)
        )
        pos_mat = permutedims(
            pyconvert(Matrix{Float64}, py_pos_mat),
            (2, 1)
        )
        pos_mat = mod1.(pos_mat, 1)
        num_vec = pyconvert(Vector{Int64}, py_num_vec)

        mag_flag_vec = [(num in mag_num_vec) for num in num_vec]
        mag_count = sum(mag_flag_vec)
        mag_config_array = mag_config(mag_count, target_idx_vec)

        struc_vec = Struc[]
        for idx in axes(mag_config_array, 3)
            spin_mat = zero(pos_mat)
            spin_mat[:, mag_flag_vec] .= mag_config_array[:, :, idx]

            struc = Struc(
                idx,
                lattice_mat,
                num_vec,
                pos_mat,
                spin_mat
            )

            push!(struc_vec, struc)
        end

        return struc_vec
    end


    function equal_pair(py_struc, spg_num, mag_num_vec, target_idx_vec)
        # remove all nonmagnetic atoms for only considering pairs between
        # magnetic atoms
        py_mag_struc = py_struc.copy()
        num_vec = pyconvert(Vector, py_struc.atomic_numbers)
        nonmag_idx_vec = findall([!(num in mag_num_vec) for num in num_vec]) .- 1
        py_mag_struc.remove_sites(PyList(nonmag_idx_vec))

        cutoff = ceil(pyconvert(
            Float64,
            py_mag_struc.get_distance((target_idx_vec .- 1)...)
        ))  # get the length of given atom pair
        (
            py_center_indices,
            py_points_indices,
            _,
            _,
            py_symmetry_indices,
            _
        ) = py_mag_struc.get_symmetric_neighbor_list(
            cutoff,
            spg_num,
            numerical_tol=1e-2,
            exclude_self=true,
            unique=false
        )
        center_idx_vec = pyconvert(Vector, py_center_indices) .+ 1
        points_idx_vec = pyconvert(Vector, py_points_indices) .+ 1
        symmetry_idx_vec = pyconvert(Vector{Int64}, py_symmetry_indices)

        target_center_idx_vec = findall(==(target_idx_vec[1]), center_idx_vec)
        target_points_idx_vec = findall(==(target_idx_vec[2]), points_idx_vec)
        target_pair_idx_vec = intersect(target_center_idx_vec, target_points_idx_vec)
        if length(target_pair_idx_vec) > 1
            error("Given supercell is not large enough for target atom pair.")
        else
            target_pair_idx = target_pair_idx_vec[1]
        end

        sym_idx = symmetry_idx_vec[target_pair_idx]
        equal_pair_idx_vec = intersect(
            target_center_idx_vec,
            findall(==(sym_idx), symmetry_idx_vec)
        )
    end


    function to_vasp_inputs(
        map::Map;
        incar_path="./INCAR",
        poscar_path="./POSCAR",
        potcar_path="./POTCAR",
        kpoints_path="./KPOINTS"
    )
        par_dir_name = "./J_MAT/"
        if isdir(par_dir_name)
            @warn "Old input dir detected! Do you want to delete them? (Y[es]/N[o])"
            choice = readline(stdin)
            if uppercase(first(choice)) == 'Y'
                rm(par_dir_name, recursive=true, force=true)
            else
                @info "Do not delete them."
                return nothing
            end
        end

        py_incar = py_Incar.from_file(incar_path)

        mkdir(par_dir_name)
        conf_dir_vec = String[]
        for (struc_idx, struc) in enumerate(map.struc_vec)
            conf_dir = par_dir_name * "conf_$(struc_idx)/"
            mkdir(conf_dir)

            py_magmom_list = py_np.array(transpose(struc.spin_mat)).tolist()
            py_incar.update(Dict(
                "MAGMOM" => py_magmom_list,
                "M_CONSTR" => pylist(struc.spin_mat)
            ))

            py_incar.write_file(conf_dir * "INCAR")
            symlink(abspath(poscar_path), conf_dir * "POSCAR")
            symlink(abspath(potcar_path), conf_dir * "POTCAR")
            symlink(abspath(kpoints_path), conf_dir * "KPOINTS")

            push!(conf_dir_vec, conf_dir)
        end

        @info "Storing path to different configuration into `CONF_DIR`. One may use SLURM's job array to calculate."
        writedlm("CONF_DIR", conf_dir_vec)
    end
end