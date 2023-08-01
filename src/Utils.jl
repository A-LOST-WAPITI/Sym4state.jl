module Utils
    using InvertedIndices
    using LinearAlgebra
    using PythonCall
    using ..Types

    
    export fourstate, simplify_map


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

    function simplify_map(fallback::FallbackList)
        map_array = [zeros(Int8, 4) for _ = 1:3, _ = 1:3]
        temp = eachslice(reshape(fallback.(1:36), (4, 3, 3)), dims=(2, 3))

        for idx in eachindex(temp)
            element_comp = temp[idx]
            part1 = element_comp[[1, 4]]
            part2 = element_comp[[2, 3]]

            part_diff = setdiff(part1, part2)
            if length(part_diff) != 0
                map_array[idx] .= element_comp
            end
        end

        return map_array
    end


end