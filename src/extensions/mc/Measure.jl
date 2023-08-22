module MCMeasure
    using Statistics
    using StatsBase

    function measure_mag_mean(states_array)
        mag_mean_vec = mean(states_array, dims=(1, 2, 3))
        mag_mean = norm(mag_mean_vec)

        return mag_mean
    end

    function site_energy!(
        energy_array,
        states_array,
        mag_field,
        a_mat,
        j_array,
        idx_mat,
        atom_type
    )
        x = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        y = (blockIdx().y - 1) * blockDim().y + threadIdx().y

        n_x, n_y, _, _ = size(states_array)
        neigh_type_num = size(j_array, 3)
        if x <= n_x && y <= n_y
            state_vec = @view states_array[x, y, atom_type, :]

            energy::Float32 = tri_product(state_vec, a_mat, state_vec) - (
                mag_field[1] * state_vec[1] +
                mag_field[2] * state_vec[2] +
                mag_field[3] * state_vec[3]
            )
            for i = 1:neigh_type_num
                j_mat = @view j_array[:, :, i]
                x_diff, y_diff, target_type = @view idx_mat[i, :]
                x_target, y_target = target_idx(
                    x,
                    y,
                    x_diff,
                    y_diff,
                    n_x,
                    n_y
                )
                neigh_state_vec = @view states_array[x_target, y_target, target_type, :]

                energy += tri_product(
                    state_vec,
                    j_mat,
                    neigh_state_vec
                )
            end

            energy_array[x, y, atom_type] = energy
        end

        return nothing
    end

    function measure_energy(
        states_array,
        mag_field_tuple,
        a_tuple,
        j_tuple,
        idx_tuple,
        threads_tuple,
        blocks_tuple
    )
        n_x, n_y, atom_type_num, _ = size(states_array)
        energy_array = CUDA.zeros(n_x, n_y, atom_type_num)
        for (
            atom_type, 
            (mag_field, a_mat, j_array, idx_mat)
        ) in enumerate(zip(mag_field_tuple, a_tuple, j_tuple, idx_tuple))
            @cuda threads=threads_tuple blocks=blocks_tuple site_energy!(
                energy_array,
                states_array,
                mag_field,
                a_mat,
                j_array,
                idx_mat,
                atom_type
            )
        end

        return sum(energy_array)/2
    end
end