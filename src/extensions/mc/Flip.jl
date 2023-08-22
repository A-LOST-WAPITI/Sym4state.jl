module MCFlip
    using LinearAlgebra
    using CUDA
    using PhysicalConstants.CODATA2018  # 参考CODATA2018标准的物理常数
    using Unitful           # 单位运算软件包
    using UnitfulAtomic     # 原子单位运算软件包
    using Markdown


    const KB::Float32 = BoltzmannConstant |> auconvert |> austrip   # 玻尔兹曼常数


    function target_idx(x, y, x_diff, y_diff, n_x, n_y)
        x_target = mod1(x + x_diff, n_x)
        y_target = mod1(y + y_diff, n_y)

        return x_target, y_target
    end

    function tri_product(s_vec, j_mat, t_vec)
        result = 0.0f0

        for i = 1:3
            result += t_vec[i] * (
                j_mat[1, i] * s_vec[1] +
                j_mat[2, i] * s_vec[2] +
                j_mat[3, i] * s_vec[3]
            )
        end

        return result
    end

    function fixed_rand_states(n_x, n_y, atom_type_num)
        rand_states_array = CUDA.zeros(n_x, n_y, atom_type_num, 3)
        z2_array = CUDA.rand(n_x, n_y, atom_type_num)
        phi_array = CUDA.rand(n_x, n_y, atom_type_num) * pi

        qx = @. sqrt(z2_array) * sin(phi_array)
        qy = @. sqrt(z2_array) * cos(phi_array)
        qz = @. sqrt(1 - z2_array) * sin(phi_array)
        qw = @. sqrt(1 - z2_array) * cos(phi_array)

        rand_states_array[:, :, :, 1] .= @. 2qy * qw - 2qx * qz
        rand_states_array[:, :, :, 2] .= @. 2qy * qz + 2qx * qw
        rand_states_array[:, :, :, 3] .= @. 1 - 2qx^2 - 2qy^2

        return rand_states_array
    end

    function try_flip_1d!(
        states_array,
        rand_states_array,
        prob_array,
        temperature,
        mag_field,
        a_mat,
        j_array,
        idx_mat,
        atom_type,
        para_dis,
        choice
    )
        x = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        y = (blockIdx().y - 1) * blockDim().y + threadIdx().y

        n_x, n_y, _, _ = size(states_array)
        neigh_type_num = size(j_array, 3)
        if x <= n_x && y <= n_y && (x + y)%para_dis == choice
            state_vec = @view states_array[x, y, atom_type, :]
            rand_state_vec = @view rand_states_array[x, y, atom_type, :]

            energy_raw::Float32 = tri_product(state_vec, a_mat, state_vec) - (
                mag_field[1] * state_vec[1] +
                mag_field[2] * state_vec[2] +
                mag_field[3] * state_vec[3]
            )
            energy_try::Float32 = tri_product(rand_state_vec, a_mat, rand_state_vec) - (
                mag_field[1] * rand_state_vec[1] +
                mag_field[2] * rand_state_vec[2] +
                mag_field[3] * rand_state_vec[3]
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

                energy_raw += tri_product(
                    state_vec,
                    j_mat,
                    neigh_state_vec
                )
                energy_try += tri_product(
                    rand_state_vec,
                    j_mat,
                    neigh_state_vec
                )
            end

            delta_energy = energy_try - energy_raw  # 计算反转产生的能量差值
            if delta_energy < 0 || prob_array[x, y, atom_type] < exp(-delta_energy/(KB * temperature))
                # 不要觉得这个操作很傻，单一CUDA kernel只支持标量运算
                for pos = 1:3
                    states_array[x, y, atom_type, pos] = rand_states_array[x, y, atom_type, pos]
                end
            end
        end

        return nothing
    end

    function try_flip_2d!(
        states_array,
        rand_states_array,
        prob_array,
        temperature,
        mag_field,
        a_mat,
        j_array,
        idx_mat,
        atom_type,
        x_win,
        y_win,
        x_choice,
        y_choice
    )
        x = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        y = (blockIdx().y - 1) * blockDim().y + threadIdx().y

        n_x, n_y, _, _ = size(states_array)
        neigh_type_num = size(j_array, 3)
        if x <= n_x && y <= n_y && x%x_win == x_choice && y%y_win == y_choice
            state_vec = @view states_array[x, y, atom_type, :]
            rand_state_vec = @view rand_states_array[x, y, atom_type, :]

            energy_raw::Float32 = tri_product(state_vec, a_mat, state_vec) - (
                mag_field[1] * state_vec[1] +
                mag_field[2] * state_vec[2] +
                mag_field[3] * state_vec[3]
            )
            energy_try::Float32 = tri_product(rand_state_vec, a_mat, rand_state_vec) - (
                mag_field[1] * rand_state_vec[1] +
                mag_field[2] * rand_state_vec[2] +
                mag_field[3] * rand_state_vec[3]
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

                energy_raw += tri_product(
                    state_vec,
                    j_mat,
                    neigh_state_vec
                )
                energy_try += tri_product(
                    rand_state_vec,
                    j_mat,
                    neigh_state_vec
                )
            end

            delta_energy = energy_try - energy_raw  # 计算反转产生的能量差值
            if delta_energy < 0 || prob_array[x, y, atom_type] < exp(-delta_energy/(KB * temperature))
                # 不要觉得这个操作很傻，单一CUDA kernel只支持标量运算
                for pos = 1:3
                    states_array[x, y, atom_type, pos] = rand_states_array[x, y, atom_type, pos]
                end
            end
        end

        return nothing
    end

    function mcmc!(
        states_array,
        temperature,
        mag_field_tuple,
        a_tuple,
        j_tuple,
        idx_tuple,
        reducible_tuple,
        threads_tuple,
        blocks_tuple
    )
        n_x, n_y, atom_type_num, _ = size(states_array)
        rand_states_array = fixed_rand_states(n_x, n_y, atom_type_num)
        prob_array = CUDA.rand(n_x, n_y, atom_type_num)

        for (
            atom_type, 
            (mag_field, a_mat, j_array, idx_mat, reducible)
        ) in enumerate(zip(mag_field_tuple, a_tuple, j_tuple, idx_tuple, reducible_tuple))
            para_dis, x_win, y_win = reducible

            if iszero(para_dis)
                for x_choice = 0:(x_win - 1), y_choice = 0:(y_win - 1)
                    @cuda threads=threads_tuple blocks=blocks_tuple try_flip_2d!(
                        states_array,
                        rand_states_array,
                        prob_array,
                        temperature,
                        mag_field,
                        a_mat,
                        j_array,
                        idx_mat,
                        atom_type,
                        x_win,
                        y_win,
                        x_choice,
                        y_choice
                    )
                end
            else
                for choice = 0:(para_dis - 1)
                    @cuda threads=threads_tuple blocks=blocks_tuple try_flip_1d!(
                        states_array,
                        rand_states_array,
                        prob_array,
                        temperature,
                        mag_field,
                        a_mat,
                        j_array,
                        idx_mat,
                        atom_type,
                        para_dis,
                        choice
                    )
                end
            end
        end

        return nothing
    end
end