module Utils
    using InvertedIndices
    using LinearAlgebra
    using DelimitedFiles
    using PythonCall
    using Printf
    using CellListMap: neighborlist
    using DataStructures: DisjointSets, find_root!, num_groups
    using ..Types
    using ..Python


    export get_all_j_struc_vec, to_vasp_inputs, equal_pair, get_py_struc, get_sym_op_vec, get_j_mat, struc_compare
    export linear_idx_to_vec
    export py_struc_to_struc


    include("data/CovalentRadius.jl")


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


    function py_struc_to_struc(py_struc; idx::Int=1, spin_mat::Union{Nothing, AbstractMatrix}=nothing)
        py_lattice_mat = py_struc.lattice.matrix.transpose()
        py_pos_mat = py_struc.frac_coords.transpose()
        py_num_vec = py_struc.atomic_numbers

        lattice_mat = pyconvert(Matrix{Float64}, py_lattice_mat)
        pos_mat = pyconvert(Matrix{Float64}, py_pos_mat)
        pos_mat = mod1.(pos_mat, 1)
        num_vec = pyconvert(Vector{Int64}, py_num_vec)
        if isa(spin_mat, Nothing)
            spin_mat = zero(pos_mat)
        end

        struc = Struc(
            idx,
            lattice_mat,
            num_vec,
            pos_mat,
            spin_mat
        )

        return struc
    end


    function get_all_j_struc_vec(struc::Struc, mag_num_vec, target_idx_vec)
        num_vec = struc.num_vec

        mag_flag_vec = [(num in mag_num_vec) for num in num_vec]
        mag_count = sum(mag_flag_vec)
        mag_config_array = mag_config(mag_count, target_idx_vec)

        struc_vec = Struc[]
        for idx in axes(mag_config_array, 3)
            spin_mat = zeros(3, length(num_vec))
            spin_mat[:, mag_flag_vec] .= mag_config_array[:, :, idx]

            struc = Struc(
                idx,
                struc.lattice_mat,
                struc.num_vec,
                struc.pos_mat,
                spin_mat
            )

            push!(struc_vec, struc)
        end

        return struc_vec
    end

    # TODO: generate structures used in A matrix
    function get_all_a_struc_vec_nosym(py_struc, mag_num_vec, target_idx_vec)
        
    end

    function struc_compare(x::Struc, y::Struc; atol=1e-2)
        approx_flag = true

        corresponding_dict = Dict{Int, Int}()
        for (one_atom_idx, one_atom) in enumerate(x)
            occur_flag = false
            for (another_atom_idx, another_atom) in enumerate(y)
                if isapprox(one_atom, another_atom, atol=atol)
                    occur_flag = true
                    corresponding_dict[another_atom_idx] = one_atom_idx
                    break
                end
            end

            if !occur_flag
                approx_flag = false
                break
            end
        end

        return approx_flag, corresponding_dict
    end

    function consider_idx_in_radius(struc::Struc, center_idx, cutoff_radius)
        frac_pos_mat = struc.pos_mat
        cell_mat = struc.lattice_mat
        cart_pos_mat = cell_mat * frac_pos_mat
        neighbor_list = neighborlist(
            eachcol(cart_pos_mat),
            cutoff_radius,
            unitcell=cell_mat,
            parallel=true
        )

        consider_idx_vec = Int[]
        for pair in neighbor_list
            bond_ends_set = Set(pair[1:2])
            if center_idx in bond_ends_set
                push!(
                    consider_idx_vec,
                    pop!(setdiff(bond_ends_set, center_idx))
                )
            end
        end

        return consider_idx_vec
    end

    function linear_idx_to_vec(idx, supercell_size, num_pri_sites)
        temp = reshape(1:9, num_pri_sites, reverse(supercell_size)...)

        idx_vec = findfirst(==(idx), temp).I |> collect

        reverse!(idx_vec)
        idx_vec .-= 1
        idx_vec[end] += 1   # there is no 0th atom

        for idx_pos = 1:3
            if idx_vec[idx_pos] > supercell_size[idx_pos]/2
                idx_vec[idx_pos] -= supercell_size[idx_pos]
            end
        end

        return idx_vec
    end

    function equal_pair(
        py_struc,
        mag_num_vec,
        center_idx,
        cutoff_radius,
        sym_op_vec;
        atol=1e-2
    )
        # remove all nonmagnetic atoms for only considering pairs between
        # magnetic atoms
        py_mag_struc = py_struc.copy()
        num_vec = pyconvert(Vector, py_struc.atomic_numbers)
        nonmag_idx_vec = setdiff(unique(num_vec), mag_num_vec)
        py_mag_struc.remove_species(PyList(nonmag_idx_vec))

        mag_struc = py_struc_to_struc(py_mag_struc)
        consider_idx_vec = consider_idx_in_radius(mag_struc, center_idx, cutoff_radius)

        # there exists multiple pairs between center atom and one point atom
        if length(unique(consider_idx_vec)) != length(consider_idx_vec)
            error(
                "Given supercell is not large enough " *
                "for calculating all interactions within a cutoff radius of $(cutoff_radius) Å.")
        end

        pair_ds = DisjointSets(consider_idx_vec)
        pair_relation_dict = Dict{AbstractVector{Int}, SymOp}()
        for sym_op in sym_op_vec
            mag_struc_after_op = sym_op * mag_struc

            # without considering spin, those two structures should be approximate
            _, corresponding_dict = struc_compare(
                mag_struc,
                mag_struc_after_op,
                atol=atol
            )

            center_after_op_idx = corresponding_dict[center_idx]
            if center_after_op_idx != center_idx    # fix the center point
                continue
            end
            for idx in consider_idx_vec # find if there is any pair
                raw_idx = corresponding_dict[idx]
                if raw_idx in consider_idx_vec
                    union!(
                        pair_ds,
                        raw_idx,
                        idx
                    )
                    pair_vec = [raw_idx, idx]
                    if !haskey(pair_relation_dict, pair_vec)
                        pair_relation_dict[pair_vec] = sym_op
                    end
                end
            end
        end

        return pair_ds, pair_relation_dict
    end


    function set_rwigs(poscar_path)
        symbol_vec = split(readlines(poscar_path)[6])

        rwigs_vec = Float64[]
        for element_symbol in symbol_vec
            item = RADIUS_DICT[element_symbol]
            names = item["name"]
            radius = item["radii"]

            mult = length(names)
            if mult == 1
                radii = radius[1]
            else
                @info "There are different hybridisations of this element."
                @info "     Hybridisation: " * repr(names)
                @info "Coorespoding radii: " * repr(radius)
                @info "Type the index of the radii you want:"
                idx = parse(Int64, readline(stdin))
                
                radii = radius[idx]
            end

            @info "RWIGS of $(element_symbol) is set to $(radii)."
            push!(rwigs_vec, radii)
        end

        return rwigs_vec
    end


    function to_vasp_inputs(
        center_map_vec::Vector{Vector{Map}};
        incar_path="./INCAR",
        poscar_path="./POSCAR",
        potcar_path="./POTCAR",
        kpoints_path="./KPOINTS",
        kwargs...
    )
        @assert isfile(incar_path) "Invalid `INCAR` path!"
        @assert isfile(potcar_path) "Invalid `POTCAR` path!"
        @assert isfile(kpoints_path) "Invalid `KPOINTS` path!"

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
        rwigs_vec = set_rwigs(poscar_path)

        mkpath(par_dir_name)
        conf_dir_vec = String[]
        for (center_idx, map_vec) in enumerate(center_map_vec)
            center_dir_name = par_dir_name * "center_$(center_idx)/"
            for (group_idx, map) in enumerate(map_vec)
                group_dir_name = center_dir_name * "goup_$(group_idx)/"
                for (struc_idx, struc) in enumerate(map.struc_vec)
                    conf_dir = group_dir_name * "conf_$(struc_idx)/"
                    mkpath(conf_dir)

                    py_magmom_list = py_np.array(transpose(struc.spin_mat)).tolist()
                    update_pack_dict = Dict(
                        "MAGMOM" => py_magmom_list,
                        "M_CONSTR" => pylist(struc.spin_mat),
                        "I_CONSTRAINED_M" => 1,
                        "RWIGS" => pylist(rwigs_vec)
                    )
                    # other parameters
                    for pair in kwargs
                        push!(update_pack_dict, pair)
                    end
                    py_incar.update(update_pack_dict)

                    py_incar.write_file(conf_dir * "INCAR")
                    # using relative path for flexibility
                    symlink(relpath(poscar_path, conf_dir), conf_dir * "POSCAR")
                    symlink(relpath(potcar_path, conf_dir), conf_dir * "POTCAR")
                    symlink(relpath(kpoints_path, conf_dir), conf_dir * "KPOINTS")

                    push!(conf_dir_vec, conf_dir)
                end
            end
        end

        @info "Storing path to different configuration into `J_CONF_DIR`. One may use SLURM's job array to calculate."
        writedlm("J_CONF_DIR", conf_dir_vec)
    end


    get_py_struc(filepath::String) = py_Struc.from_file(filepath)


    function get_sym_op_vec(
        py_struc,
        mag_num_vec,
        supercell_size;
        symprec=1e-2,
        angle_tolerance=5.0
    )
        py_sga = py_Sga(
            py_struc,
            symprec=symprec,
            angle_tolerance=angle_tolerance
        )
        py_refined_struc = py_sga.get_refined_structure()
        py_num_vec = py_refined_struc.atomic_numbers
        num_vec = pyconvert(Vector{Int64}, py_num_vec)
        mag_atom_count = count(x -> x ∈ mag_num_vec, num_vec)

        py_refined_struc.make_supercell(supercell_size)
        py_sga = py_Sga(
            py_refined_struc,
            symprec=symprec,
            angle_tolerance=angle_tolerance
        )
        
        py_sym_dict = py_sga.get_symmetry_dataset()

        spg_num = pyconvert(Int64, py_sym_dict["number"])
        lattice_mat = permutedims(
            pyconvert(Array{Float64}, py_refined_struc.lattice.matrix),
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

        return spg_num, mag_atom_count, op_vec, py_refined_struc
    end


    function get_energy_and_magmom(py_outcar)
        energy = pyconvert(Float64, py_outcar.final_energy)

        py_mag_mat = py_np.vstack(
            [item["tot"].moment for item in py_outcar.magnetization]
        )
        mag_mat = permutedims(
            pyconvert(Matrix{Float64}, py_mag_mat),
            (2, 1)
        )

        return energy, mag_mat
    end

    function get_j_mat(
        map::Map,
        conf_list_path::String,
        target_idx_vec::Vector{Int64}
    )
        fallback_vec = map.fallback_vec
        conf_dir_vec = readdlm(conf_list_path)
        @assert length(fallback_vec) == length(conf_dir_vec)

        energy_vec = Dict{Int8, Float64}()
        for (conf_idx, conf_dir) in enumerate(conf_dir_vec)
            target_conf_mag_mat = map.struc_vec[conf_idx].spin_mat[:, target_idx_vec]
            py_outcar = py_Outcar(conf_dir * "OUTCAR")
            energy, mag_mat = get_energy_and_magmom(py_outcar)

            target_mag_mat = mag_mat[:, target_idx_vec]
            target_mag_mat = round.(target_mag_mat, RoundToZero; digits=1)
            for col in eachcol(target_mag_mat)
                normalize!(col)
            end
            if !isapprox(target_mag_mat, target_conf_mag_mat)
                @error "MAGMOM of $(conf_dir) changed after SCF!"
            end

            energy_vec[fallback_vec[conf_idx]] = energy
        end

        j_mat = zeros(3, 3)
        for i = 1:3, j = 1:3
            map_vec = map.map_mat[i, j]
            
            if map_vec[1] == 0
                j_mat[i, j] = 0
            else
                idx_1, idx_2, idx_3, idx_4 = map_vec
                j_mat[i, j] = (
                    energy_vec[idx_1] -
                    energy_vec[idx_2] -
                    energy_vec[idx_3] +
                    energy_vec[idx_4]
                )/4
            end
        end

        return j_mat
    end
end