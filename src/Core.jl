module ModCore
    using ProgressMeter
    using ..Pymatgen
    using ..Spglib
    using ..Utils
    using ..Types


    export sym4state


    function sym4state(
        filepath,
        mag_num_vec,
        target_idx_vec;
        atol=1e-2,
        symprec=1e-2
    )
        @info "Absolute tolrance is set to $(atol)"
        @info "Symmetry precision is set to $(symprec)"
        py_struc = get_py_struc(filepath)

        #TODO: Check whether `target_idx_vec` is enough.

        spg_num, sym_op_vec = get_sym_op_vec(py_struc, symprec=symprec)
        @info "The space group number of given structure is $(spg_num) with given `symprec`"
        struc_vec = fourstate(py_struc, mag_num_vec, target_idx_vec)

        unique_struc_vec = Struc[]
        fallback = FallbackList(36)
        p = Progress(length(struc_vec) * length(sym_op_vec))
        for struc in struc_vec
            for sym_op in sym_op_vec
                struc_after_op = sym_op * struc
                source_uni_num = struc_after_op.uni_num

                occur_flag = false
                for struc_occur in unique_struc_vec
                    target_uni_num = struc_occur.uni_num
                    if source_uni_num != target_uni_num && isapprox(
                        struc_after_op,
                        struc_occur,
                        atol=atol
                    )
                        occur_flag = true
                        fallback(
                            source_uni_num,
                            target_uni_num
                        )
                    end
                end

                if !occur_flag
                    push!(
                        unique_struc_vec,
                        struc_after_op
                    )
                end

                next!(p)
            end
        end

        return fallback
    end
end