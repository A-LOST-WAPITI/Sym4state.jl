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
        py_refined_struc.to("POSCAR")
        @info "The refined structure has been dumped into \"POSCAR\""
        
        center_idx_vec = [
            (idx - 1) * prod(supercell_size) + idx
            for idx = 1:mag_atom_count
        ]
        pair_ds, pair_relation_dict = equal_pair(
            py_refined_struc,
            mag_num_vec,
            center_idx_vec,
            cutoff_radius,
            sym_op_vec
        )

        ngroups = num_groups(pair_ds)
        @info "There are $(ngroups) different type(s) of pairs."
        pair_vec_vec = pair_ds.revmap
        parents_vec = [
            find_root!(pair_ds, pair_vec)
            for pair_vec in pair_vec_vec
        ]
        group_parents_vec = unique(parents_vec)
        for (parent_idx, group_parent) in enumerate(group_parents_vec)
            @info "For the $(parent_idx)th group, equivalent pairs are shown as follows:"
            for (pair_vec, parent) in zip(pair_vec_vec, parents_vec)
                if parent == group_parent
                    @info "$(parent) <=> $(pair_vec)"
                end
            end
        end

        raw_struc = py_struc_to_struc(py_refined_struc)
        map_vec = Map[]
        relation_vec = CoeffMatRef[]
        map::Union{Nothing, Map} = nothing
        for (group_idx, group_parent) in enumerate(group_parents_vec)
            min_energy_num = 37 # make sure `min_energy_num` can be updated
            min_pair_vec = Int[]
            for (pair_vec, parent) in zip(pair_vec_vec, parents_vec)
                # pairs from different groups are skipped
                if parent != group_parent
                    continue
                end

                println()
                @info "Checking pair $(pair_vec) ..."
                temp_struc_vec = get_all_j_struc_vec(raw_struc, mag_num_vec, pair_vec)
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

                    min_pair_vec = pair_vec
                    temp_map = Map(
                        temp_fallback_ds,
                        temp_struc_vec
                    )
                    min_energy_num = length(temp_map.fallback_vec)
                    map = temp_map

                    @info "The number of energies now is $(min_energy_num)."
                end
            end

            push!(map_vec, map)
            for (pair_vec, parent) in zip(pair_vec_vec, parents_vec)
                if parent != group_parent
                    continue
                end

                pair_relation_vec = vcat(
                    min_pair_vec,
                    pair_vec
                )
                center_idx = linear_idx_to_vec(pair_vec[1], supercell_size, mag_atom_count)
                point_idx = linear_idx_to_vec(pair_vec[2], supercell_size, mag_atom_count)

                cell_idx_diff = point_idx[2:end] .- center_idx[2:end]
                fixed_pair_vec = vcat(center_idx[1], cell_idx_diff[1:2], point_idx[1])
                coeff_ref = CoeffMatRef(
                    group_idx,
                    fixed_pair_vec,
                    pair_relation_dict[pair_relation_vec]
                )
                push!(
                    relation_vec,
                    coeff_ref
                )
            end
        end

        return map_vec, relation_vec
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
        map_vec, relation_vec = sym4state(
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
        @info "Saving the reduced map and relations into \"cal.jld2\"..."
        save(
            "cal.jld2",
            Dict(
                "map" => map_vec,
                "relation" => relation_vec
            )
        )

        to_vasp_inputs(
            map_vec,
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
        map_vec = load(
            filepath,
            "map"
        )

        to_vasp_inputs(
            map_vec,
            incar_path=incar_path,
            poscar_path=poscar_path,
            potcar_path=potcar_path,
            kpoints_path=kpoints_path,
            kwargs...
        )
    end

    function post_process(
        cal_file_path::String
    )
        cal_dir = dirname(abspath(cal_file_path)) * "/cal/"

        map_vec, relation_vec = load(
            cal_file_path,
            "map",
            "relation"
        )

        pair_mat, coeff_array = get_pair_and_coeff(
            map_vec,
            relation_vec,
            cal_dir
        )

        return pair_mat, coeff_array
    end
end