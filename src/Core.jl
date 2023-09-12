module ModCore
    using ProgressMeter
    using FileIO
    using DataStructures: IntDisjointSets, union!, num_groups, find_root!
    using ..Utils
    using ..Types


    export pre_process


    function sym4state(
        py_struc,
        mag_num_vec,
        supercell_size,
        cutoff_radius;
        atol=1e-2,
        symprec=1e-2,
        angle_tolerance=5.0,
        show_progress_bar=true,
        center_idx=1
    )
        @info "Absolute tolrance is set to $(atol)"
        @info "Symmetry precision is set to $(symprec)"

        spg_num, sym_op_vec, py_refined_struc = get_sym_op_vec(
            py_struc,
            supercell_size,
            symprec=symprec,
            angle_tolerance=angle_tolerance
        )
        @info "The space group number of given structure is $(spg_num) with given `symprec`"
        dump_struc_name = "POSCAR_refined_$(supercell_size[1])x$(supercell_size[2])"    # we only consider 2D materials
        py_refined_struc.to(dump_struc_name)
        @info "The refined structure has been dumped into \"$(dump_struc_name)\""
        
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

        min_energy_num = 37 # make sure `min_energy_num` can be updated
        min_point_idx = 0
        struc_vec = Struc[]
        fallback_ds = IntDisjointSets(37)   # make sure `fallback_ds` can be updated
        for point_idx in pair_ds.revmap
            @info "Checking pair $(center_idx) <=> $(point_idx) ..."
            target_idx_vec = [center_idx, point_idx]
            temp_struc_vec = get_all_j_struc_vec(py_refined_struc, mag_num_vec, target_idx_vec)
            mag_struc_vec = [magonly(struc, mag_num_vec) for struc in temp_struc_vec]

            @info "Reducing 4-state J matrix ..."
            unique_mag_struc_vec = Struc[]
            temp_fallback_ds = IntDisjointSets(36)
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
                        approx_flag, _ = struc_compare(
                            mag_struc_after_op,
                            mag_struc_occur,
                            atol=atol
                        )
                        if source_uni_num != target_uni_num && approx_flag
                            Threads.lock(para_lock) do
                                occur_flag = true

                                union!(
                                    temp_fallback_ds,
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

            energy_num = num_groups(temp_fallback_ds)
            if energy_num < min_energy_num
                @info "Find pair with higher symmetry!"
                @info "The number of energies now is $(energy_num)."

                min_energy_num = energy_num
                min_point_idx = point_idx
                struc_vec = temp_struc_vec
                fallback_ds = temp_fallback_ds
            end
        end

        map = Map(fallback_ds, struc_vec)

        @info "Saving the reduced map into \"Map.jld2\"..."
        save(
            "Map.jld2",
            Dict(
                "map" => map
            )
        )

        return map
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
        center_idx=1,
        incar_path="./INCAR",
        potcar_path="./POTCAR",
        kpoints_path="./KPOINTS"
    )
        py_struc = get_py_struc(filepath)
        map = sym4state(
            py_struc,
            mag_num_vec,
            supercell_size,
            cutoff_radius;
            atol=atol,
            symprec=symprec,
            angle_tolerance=angle_tolerance,
            show_progress_bar=show_progress_bar,
            center_idx=center_idx
        )

        to_vasp_inputs(
            map,
            incar_path=incar_path,
            poscar_path=poscar_path,
            potcar_path=potcar_path,
            kpoints_path=kpoints_path
        )
    end

    function post_process(
        map_path::String,
        target_idx_vec
    )
        map = load(
            map_path,
            "map"
        )
        cal_dir = dirname(map_path)
        conf_list_path = cal_dir * "/J_CONF_DIR"

        dir_now = pwd()
        cd(cal_dir)
        j_mat = get_j_mat(map, conf_list_path, target_idx_vec)
        cd(dir_now)

        return j_mat
    end
end