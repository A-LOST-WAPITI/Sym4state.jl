module ModCore
    using ProgressMeter
    using FileIO
    using DataStructures: IntDisjointSets, union!, num_groups, find_root!
    using ..Utils
    using ..Types


    export pre_process


    function reduce_j_mat_for_a_pair(
        mag_struc_vec,
        sym_op_vec;
        show_progress_bar=true,
        atol=1e-2
    )
        unique_mag_struc_vec = Struc[]
        fallback_ds = IntDisjointSets(36)
        p = Progress(
            length(mag_struc_vec) * length(sym_op_vec);
            showspeed=true,
            enabled=show_progress_bar
        )
        for mag_struc in mag_struc_vec
            for sym_op in sym_op_vec
                mag_struc_after_op = sym_op * mag_struc
                source_uni_num = mag_struc_after_op.uni_num

                occur_flag = false
                para_lock = Threads.SpinLock()
                Threads.@threads for mag_struc_occur in unique_mag_struc_vec
                    target_uni_num = mag_struc_occur.uni_num
                    approx_flag, _ = struc_compare( # corresponding_dict not needed here
                        mag_struc_after_op,
                        mag_struc_occur,
                        atol=atol
                    )
                    if source_uni_num != target_uni_num && approx_flag
                        Threads.lock(para_lock) do
                            occur_flag = true

                            union!(
                                fallback_ds,
                                source_uni_num,
                                target_uni_num
                            )
                        end
                    end
                end

                if !occur_flag
                    push!(
                        unique_mag_struc_vec,
                        mag_struc_after_op
                    )
                end

                next!(p)
            end
        end

        return fallback_ds
    end


    function sym4state(
        py_struc,
        mag_num_vec,
        supercell_size,
        cutoff_radius;
        atol=1e-2,
        symprec=1e-2,
        angle_tolerance=5.0,
        show_progress_bar=true
    )
        @info "Absolute tolrance is set to $(atol)"
        @info "Symmetry precision is set to $(symprec)"

        spg_num, mag_atom_count, sym_op_vec, py_refined_struc = get_sym_op_vec(
            py_struc,
            mag_num_vec,
            supercell_size,
            symprec=symprec,
            angle_tolerance=angle_tolerance
        )
        @info "There are $(mag_atom_count) atoms taken as magnetic in the given primitive cell."
        @info "The space group number of given structure is $(spg_num) with given `symprec`"
        dump_struc_name = "POSCAR_refined_$(supercell_size[1])x$(supercell_size[2])x$(supercell_size[3])"    # we only consider 2D materials
        py_refined_struc.to(dump_struc_name)
        @info "The refined structure has been dumped into \"$(dump_struc_name)\""
        
        center_idx_vec = [
            (idx - 1) * prod(supercell_size) + idx
            for idx = 1:mag_atom_count
        ]
        center_map_vec = Vector{Map}[]
        for center_idx in center_idx_vec
            println()
            @info "Choose the following atom as center:"
            @info "---"
            println(py_refined_struc[center_idx - 1])   # TODO: Not so pretty ...
            @info "---"
            pair_ds, pair_relation_dict = equal_pair(
                py_refined_struc,
                mag_num_vec,
                center_idx,
                cutoff_radius,
                sym_op_vec
            )

            ngroups = num_groups(pair_ds)
            @info "There are $(ngroups) different type(s) of pairs."
            point_idx_vec = pair_ds.revmap
            pair_parents_vec = [
                find_root!(pair_ds, point_idx)
                for point_idx in point_idx_vec
            ]
            group_parents_vec = unique(pair_parents_vec)
            for (parent_idx, group_parent) in enumerate(group_parents_vec)
                @info "For the $(parent_idx)th group, equivalent pairs are shown as follows:"
                for (point_idx, pair_parent) in zip(point_idx_vec, pair_parents_vec)
                    if pair_parent == group_parent
                        @info "$(center_idx) <=> $(point_idx)"
                    end
                end
            end

            raw_struc = py_struc_to_struc(py_refined_struc)
            map_vec = Map[]
            for group_parent in group_parents_vec
                min_energy_num = 37 # make sure `min_energy_num` can be updated
                min_point_idx = 0
                struc_vec = Struc[]
                fallback_ds = IntDisjointSets(37)   # make sure `fallback_ds` can be updated
                for (point_idx, pair_parent) in zip(point_idx_vec, pair_parents_vec)
                    # pairs from different groups are skipped
                    if pair_parent != group_parent
                        continue
                    end

                    @info "Checking pair $(center_idx) <=> $(point_idx) ..."
                    target_idx_vec = [center_idx, point_idx]
                    temp_struc_vec = get_all_j_struc_vec(raw_struc, mag_num_vec, target_idx_vec)
                    mag_struc_vec = [magonly(struc, mag_num_vec) for struc in temp_struc_vec]

                    @info "Reducing 4-state J matrix ..."
                    temp_fallback_ds = reduce_j_mat_for_a_pair(
                        mag_struc_vec,
                        sym_op_vec;
                        atol=atol,
                        show_progress_bar=show_progress_bar
                    )

                    energy_num = num_groups(temp_fallback_ds)
                    # record the highest symmetry for now
                    if energy_num < min_energy_num
                        @info "Find pair with higher symmetry!"
                        @info "The number of energies now is $(energy_num)."

                        min_energy_num = energy_num
                        min_point_idx = point_idx
                        struc_vec = temp_struc_vec
                        fallback_ds = temp_fallback_ds
                    end
                end

                op_dict = Dict{Vector{Int}, SymOp}()
                for point_idx in pair_ds.revmap
                    idx_diff = linear_idx_to_vec(
                        point_idx,
                        supercell_size,
                        mag_atom_count
                    )
                    op_dict[idx_diff] = pair_relation_dict[[min_point_idx, point_idx]]
                end
                map = Map(
                    fallback_ds,
                    struc_vec,
                    op_dict
                )
                push!(map_vec, map)
            end

            push!(center_map_vec, map_vec)
        end

        return center_map_vec
    end


    function pre_process(
        filepath,
        mag_num_vec,
        supercell_size,
        cutoff_radius;
        atol=1e-2,
        symprec=1e-2,
        angle_tolerance=5.0,
        show_progress_bar=true,
        incar_path="./INCAR",
        poscar_path="./POSCAR",
        potcar_path="./POTCAR",
        kpoints_path="./KPOINTS",
        kwargs...
    )
        py_struc = get_py_struc(filepath)
        center_map_vec = sym4state(
            py_struc,
            mag_num_vec,
            supercell_size,
            cutoff_radius;
            atol=atol,
            symprec=symprec,
            angle_tolerance=angle_tolerance,
            show_progress_bar=show_progress_bar
        )

        println()
        @info "Saving the reduced map into \"Map.jld2\"..."
        save(
            "Map.jld2",
            Dict(
                "map" => center_map_vec
            )
        )


        to_vasp_inputs(
            center_map_vec,
            incar_path=incar_path,
            poscar_path=poscar_path,
            potcar_path=potcar_path,
            kpoints_path=kpoints_path,
            kwargs...
        )
    end
    function pre_process(
        filepath;
        incar_path="./INCAR",
        poscar_path="./POSCAR",
        potcar_path="./POTCAR",
        kpoints_path="./KPOINTS",
        kwargs...
    )
        center_map_vec = load(
            filepath,
            "map"
        )

        to_vasp_inputs(
            center_map_vec,
            incar_path=incar_path,
            poscar_path=poscar_path,
            potcar_path=potcar_path,
            kpoints_path=kpoints_path,
            kwargs...
        )
    end

    function post_process(
        map_file_path::String;
        cal_dir::String=dirname(abspath(map_file_path)) * "/J_MAT/"
    )
        center_map_vec = load(
            map_file_path,
            "map"
        )

        point_idx_array, interact_coeff_array = get_point_and_coeff(
            center_map_vec,
            cal_dir
        )

        return point_idx_array, interact_coeff_array
    end
end