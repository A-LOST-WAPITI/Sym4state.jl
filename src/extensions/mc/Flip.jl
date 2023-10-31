module MCFlip
    using ..MCTypes

    using KernelAbstractions: @kernel, @index, get_backend


    export rand_states!, site_energy!, get_point_states!, try_flip!

    @kernel function rand_states_kernel!(states_array)
        idx = @index(Global, Cartesian)

        z2 = rand()
        phi = rand() * pi

        qx = sqrt(z2) * sin(phi)
        qy = sqrt(z2) * cos(phi)
        qz = sqrt(1 - z2) * sin(phi)
        qw = sqrt(1 - z2) * cos(phi)

        @inbounds states_array[1, idx] = 2qy * qw - 2qx * qz
        @inbounds states_array[2, idx] = 2qy * qz + 2qx * qw
        @inbounds states_array[3, idx] = 1 - 2qx^2 - 2qy^2
    end

    function rand_states!(states_array)
        array_shape = size(states_array)
        @assert length(array_shape) == 4
        @assert array_shape[1] == 3

        site_size = array_shape[2:end]
        backend = get_backend(states_array)

        kernel! = rand_states_kernel!(backend)
        kernel!(states_array, ndrange=site_size)
    end

    @kernel function try_flip_kernel!(
        states_array,
        @Const(rand_states_array),
        @Const(pair_mat),
        @Const(interact_coeff_array),
        @Const(check_array),
        @Const(magmom_vector),
        @Const(magnetic_field),
        @Const(temperature)
    )
        idx_site = @index(Global, Cartesian)
        idx_t, idx_x, idx_y = @index(Global, NTuple)
        _, _, n_x, n_y = size(states_array)
        n_p = size(interact_coeff_array, 3)

        # only give a try on sites those could be checked parallelly
        if check_array[idx_site]
            @inbounds raw_state = @view states_array[:, idx_site]
            @inbounds try_state = @view rand_states_array[:, idx_site]

            # init two energies
            raw_energy = zero(eltype(states_array))
            try_energy = zero(eltype(states_array))

            # energy from magnetic field
            for idx_pos = 1:3
                @inbounds raw_energy -= magnetic_field[idx_pos] * raw_state[idx_pos]
                @inbounds try_energy -= magnetic_field[idx_pos] * try_state[idx_pos]
            end
            raw_energy *= magmom_vector[idx_t]
            try_energy *= magmom_vector[idx_t]
            # energy from interacting
            for idx_p = 1:n_p
                @inbounds pair_vec = @view pair_mat[:, idx_p]
                if pair_vec[1] != idx_t
                    continue
                end
                target_idx_x = mod1(idx_x + pair_vec[2], n_x)  # PBC
                target_idx_y = mod1(idx_y + pair_vec[3], n_y)  # PBC
                target_idx_t = pair_vec[4]

                @inbounds interact_coeff_mat = @view interact_coeff_array[:, :, idx_p]
                @inbounds point_state = @view states_array[:, target_idx_t, target_idx_x, target_idx_y]

                for i = 1:3, j = 1:3
                    @inbounds raw_energy += raw_state[i] * interact_coeff_mat[i, j] * point_state[j]
                    @inbounds try_energy += try_state[i] * interact_coeff_mat[i, j] * point_state[j]
                end
            end

            delta_energy = try_energy - raw_energy
            if delta_energy < 0 || rand() < exp(-delta_energy/temperature)
                for idx_pos = 1:3
                    @inbounds raw_state[idx_pos] = try_state[idx_pos]
                end
            end
        end
    end

    function try_flip!(
        states_array,
        rand_states_array,
        pair_mat,
        interact_coeff_array,
        check_array,
        magmom_vector,
        magnetic_field,
        temperature
    )
        site_size = size(states_array)[2:end]
        backend = get_backend(states_array)

        kernel! = try_flip_kernel!(backend)
        kernel!(
            states_array,
            rand_states_array,
            pair_mat,
            interact_coeff_array,
            check_array,
            magmom_vector,
            magnetic_field,
            temperature,
            ndrange=site_size
        )
    end

end