module ModCore
    using FileIO
    using DataStructures: DisjointSets, IntDisjointSets, union!, num_groups, find_root!
    using LinearAlgebra: norm
    using ..Python
    using ..Utils
    using ..Types


    export pre_process


    function reduce_interact_mat_for_a_pair(
        mag_struc_vec,
        sym_op_vec;
        atol=1e-2
    )
        mag_struc_count = length(mag_struc_vec)
        @assert mag_struc_count in [20, 36]

        unique_mag_struc_dict = Dict{Int, AbstractVector{Struc}}()
        fallback_ds = IntDisjointSets(mag_struc_count)
        for mag_struc in mag_struc_vec
            target_uni_num = mag_struc.uni_num
            occur_flag = false
            for (source_uni_num, one_group_struc_vec) in unique_mag_struc_dict
                for mag_struc_occur in one_group_struc_vec
                    approx_flag, _ = struc_compare( # corresponding_dict is not needed here
                        mag_struc,
                        mag_struc_occur,
                        atol=atol
                    )

                    if approx_flag
                        union!(
                            fallback_ds,
                            source_uni_num,
                            target_uni_num
                        )
                        occur_flag = true

                        break
                    end
                end

                occur_flag && break
            end

            if !occur_flag
                new_group_struc_vec = [
                    sym_op * mag_struc
                    for sym_op in sym_op_vec
                ]

                unique_mag_struc_dict[target_uni_num] = new_group_struc_vec
            end
        end

        return fallback_ds
    end


    function supercell_check(
        py_refined_struc,
        mag_num_vec,
        mag_atom_count,
        cutoff_radius;
        symprec=1e-2,
        angle_tolerance=5.0,
        max_supercell=10
    )
        local py_refined_supercell_struc, supercell_size    # structure
        local sym_op_vec # symmetry related
        local pair_ds, pair_relation_dict   # pair
        local pass_flag
        for supercell_ratio = 2:max_supercell
            supercell_size = [supercell_ratio, supercell_ratio, 1]

            py_refined_supercell_struc = py_refined_struc.copy()
            py_refined_supercell_struc.make_supercell(supercell_size)
            py_supercell_sga = py_Sga(
                py_refined_supercell_struc,
                symprec=symprec,
                angle_tolerance=angle_tolerance
            )

            mag_struc = magonly(py_struc_to_struc(py_refined_supercell_struc), mag_num_vec)
            sym_op_vec = get_sym_op_vec(py_supercell_sga)

            center_idx_vec = [
                (idx - 1) * prod(supercell_size) + idx
                for idx = 1:mag_atom_count
            ]
            pass_flag, pair_ds, pair_relation_dict = equal_pair(
                mag_struc,
                center_idx_vec,
                cutoff_radius,
                sym_op_vec
            )

            if pass_flag
                break
            end
        end

        if !pass_flag
            error("Can't find large enough supercell!")
        end

        @info "$(supercell_size) supercell is large enough."
        py_refined_supercell_struc.to("POSCAR_refined")
        @info "The refined structure has been dumped into \"POSCAR_refined\"."
        struc = py_struc_to_struc(py_refined_supercell_struc)

        return struc, supercell_size, sym_op_vec, pair_ds, pair_relation_dict
    end


    function sym4state(
        py_struc::Py,
        mag_num_vec::AbstractVector{Int},
        cutoff_radius::T;
        atol=1e-2,
        symprec=1e-2,
        angle_tolerance=5.0,
        max_supercell=10,
        s_value=1.0
    ) where T
        py_sga = py_Sga(
            py_struc,
            symprec=symprec,
            angle_tolerance=angle_tolerance
        )
        py_refined_struc = py_sga.get_refined_structure()
        py_num_vec = py_refined_struc.atomic_numbers
        num_vec = pyconvert(Vector{Int64}, py_num_vec)
        mag_atom_count = count(∈(mag_num_vec), num_vec)
        spg_num = pyconvert(Int, py_sga.get_space_group_number())

        @info "There are $(mag_atom_count) atoms taken as magnetic in the given primitive cell."
        @info "The space group number of given structure is $(spg_num) with given `symprec`"

        raw_struc, supercell_size, sym_op_vec, pair_ds, pair_relation_dict = supercell_check(
            py_refined_struc,
            mag_num_vec,
            mag_atom_count,
            cutoff_radius;
            symprec=symprec,
            angle_tolerance=angle_tolerance,
            max_supercell=max_supercell
        )

        sym4state(
            raw_struc,
            supercell_size,
            mag_num_vec,
            sym_op_vec,
            pair_ds,
            pair_relation_dict,
            atol=atol,
            s_value=s_value
        )
    end

    function sym4state(
        raw_struc::Struc,
        supercell_size::AbstractVector{Int},
        mag_num_vec::AbstractVector{Int},
        sym_op_vec::AbstractVector{SymOp},
        pair_ds::DisjointSets,
        pair_relation_dict::Dict;
        atol=1e-2,
        s_value=1.0
    )
        rotation_sym_flag = false
        for sym_op in sym_op_vec
            if check_z_rot(sym_op)
                rotation_sym_flag = true
                break
            end
        end

        @info ""
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

        map_vec = Map[]
        relation_vec = CoeffMatRef[]
        map::Union{Nothing, Map} = nothing
        for (group_idx, group_parent) in enumerate(group_parents_vec)
            @info ""
            @info "For the $(group_idx)th group:"

            min_energy_num = 37 # make sure `min_energy_num` can be updated
            min_pair_vec = Int[]
            for (pair_vec, parent) in zip(pair_vec_vec, parents_vec)
                # pairs from different groups are skipped
                if parent != group_parent
                    continue
                end

                temp_struc_vec = get_all_interact_struc_vec(raw_struc, mag_num_vec, pair_vec, s_value)
                mag_struc_vec = [magonly(struc, mag_num_vec) for struc in temp_struc_vec]

                temp_fallback_ds = reduce_interact_mat_for_a_pair(
                    mag_struc_vec,
                    sym_op_vec;
                    atol=atol
                )

                energy_num = num_groups(temp_fallback_ds)
                # record the highest symmetry for now
                if energy_num < min_energy_num
                    @info "Find pair $(pair_vec) with higher symmetry!"

                    min_pair_vec = pair_vec
                    temp_map = Map(
                        temp_fallback_ds,
                        temp_struc_vec;
                        rotation_symmetry_flag=rotation_sym_flag
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
                fixed_pair_vec = get_fixed_pair_vec(
                    raw_struc,
                    mag_num_vec,
                    supercell_size,
                    pair_vec
                )
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
        cutoff_radius;
        atol=1e-2,
        symprec=1e-2,
        angle_tolerance=5.0,
        max_supercell=10,
        s_value=1.0,
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
            cutoff_radius;
            atol=atol,
            symprec=symprec,
            angle_tolerance=angle_tolerance,
            max_supercell=max_supercell,
            s_value=s_value
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