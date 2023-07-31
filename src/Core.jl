module ModCore
    using ..Pymatgen
    using ..Spglib
    using ..Utils
    using ..Types


    export sym4state


    function sym4state(
        filepath,
        mag_num_vec,
        target_idx_vec;
        atol=1e-2
    )
        py_struc = get_py_struc(filepath)

        #TODO: Check whether `target_idx_vec` is enough.

        spg_num, sym_op_vec = get_sym_op_vec(py_struc)
        struc_vec = fourstate(py_struc, mag_num_vec, target_idx_vec)

        unique_struc_vec = Struc[]
        fallback = FallbackList(64)
        for struc in struc_vec
            for sym_op in sym_op_vec
                struc_after_op = sym_op * struc

                for struc_occur in unique_struc_vec
                    if isapprox(
                        struc_after_op,
                        struc_occur,
                        atol=atol
                    )
                        update!(
                            fallback,
                            struc_after_op.uni_num,
                            struc_occur.uni_num
                        )
                    else
                        push!(
                            unique_struc_vec,
                            struc_after_op
                        )
                    end
                end
            end
        end

        return fallback
    end
end