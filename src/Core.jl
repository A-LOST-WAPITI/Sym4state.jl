module ModCore
    using ProgressMeter
    using FileIO
    using DataStructures: IntDisjointSets, union!
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
            symprec=symprec,
            angle_tolerance=angle_tolerance
        )
        py_refined_struc.make_supercell(supercell_size)
        @info "The space group number of given structure is $(spg_num) with given `symprec`"
        
        pair_ds, pair_relation_dict = equal_pair(
            py_refined_struc,
            mag_num_vec,
            center_idx,
            cutoff_radius,
            sym_op_vec
        )

        # TODO: Check different pair and compare them to find the minimum
        struc_vec = get_all_j_struc_vec(py_refinded_struc, mag_num_vec, target_idx_vec)
        mag_struc_vec = [magonly(struc, mag_num_vec) for struc in struc_vec]

        @info "Reducing 4-state J matrix..."
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
                    approx_flag, _ = struc_compare(
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