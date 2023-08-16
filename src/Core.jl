module ModCore
    using ProgressMeter
    using FileIO
    using ..Utils
    using ..Types


    export sym4state


    function sym4state(
        filepath,
        supercell_size,
        mag_num_vec,
        target_idx_vec;
        atol=1e-2,
        symprec=1e-2,
        angle_tolerance=5.0
    )
        @info "Absolute tolrance is set to $(atol)"
        @info "Symmetry precision is set to $(symprec)"
        py_struc = get_py_struc(filepath)

        spg_num, sym_op_vec = get_sym_op_vec(
            py_struc,
            symprec=symprec,
            angle_tolerance=angle_tolerance
        )
        @info "The space group number of given structure is $(spg_num) with given `symprec`"

        equal_pair(py_struc, supercell_size, spg_num, mag_num_vec, target_idx_vec)

        struc_vec = fourstate(py_struc, mag_num_vec, target_idx_vec)
        mag_struc_vec = [magonly(struc, mag_num_vec) for struc in struc_vec]

        @info "Reducing 4-state matrix..."
        unique_mag_struc_vec = Struc[]
        fallback = FallbackList(36)
        p = Progress(length(mag_struc_vec) * length(sym_op_vec))
        for mag_struc in mag_struc_vec
            for sym_op in sym_op_vec
                mag_struc_after_op = sym_op * mag_struc
                source_uni_num = mag_struc_after_op.uni_num

                occur_flag = false
                para_lock = Threads.SpinLock()
                Threads.@threads for mag_struc_occur in unique_mag_struc_vec
                    target_uni_num = mag_struc_occur.uni_num
                    if source_uni_num != target_uni_num && isapprox(
                        mag_struc_after_op,
                        mag_struc_occur,
                        atol=atol
                    )
                        Threads.lock(para_lock) do
                            occur_flag = true

                            fallback(
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

        map = Map(fallback, struc_vec)

        @info "Saving the reduced map into \"Map.jld2\"..."
        save(
            "Map.jld2",
            Dict(
                "map" => map
            )
        )

        return nothing
    end
end