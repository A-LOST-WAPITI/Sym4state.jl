module QMeasure
    using LinearAlgebra


    export q_value


    function q_density(si, sj, sk)
        return 2atan(
            si ⋅ (sj × sk)/(
                1 + si ⋅ sj + sj ⋅ sk + sk ⋅ si
            )
        )
    end

    function q_value(states_array::AbstractArray{T, 4}) where T
        temp_states_array = Array(states_array) # use Array
        n_x, n_y, n_t, _ = size(states_array)
        q_density_array = zeros(T, n_x, n_y, n_t)
        for idx_t = 1:n_t
            type_states_array = @view temp_states_array[idx_t, :, :, :]
            for idx_x = 1:idx_x, idx_y = 1:n_y
                idx_x_p = mod1(idx_x - 1, n_x)
                idx_x_n = mod1(idx_x + 1, n_x)
                idx_y_p = mod1(idx_y - 1, n_y)
                idx_y_n = mod1(idx_y + 1, n_y)

                q_density_array[idx_t, idx_x, idx_y] = q_density(
                    type_states_array[idx_x, idx_y, :],
                    type_states_array[idx_x_n, idx_y, :],
                    type_states_array[idx_x_n, idx_y_n, :],
                ) + q_density(
                    type_states_array[idx_x, idx_y, :],
                    type_states_array[idx_x_p, idx_y, :],
                    type_states_array[idx_x_p, idx_y_p, :],
                )
            end
        end

        q_value_vec = sum(q_density_array, dims=(2, 3))/(4 * pi)

        return q_value_vec
    end
end