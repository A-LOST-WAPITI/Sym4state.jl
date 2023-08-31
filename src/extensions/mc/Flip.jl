module MCFlip
    using ..MCTypes

    using CUDA
    using KernelAbstractions


    if CUDA.functional()
        to_gpu_or_not_to_gpu(x::AbstractArray) = CuArray(x)
    else
        to_gpu_or_not_to_gpu(x::AbstractArray) = x
    end

    @kernel function rand_states_kernel!(states_array)
        idx = @index(Global, Cartesian)

        z2 = rand()
        phi = rand() * pi

        qx = sqrt(z2) * sin(phi)
        qy = sqrt(z2) * cos(phi)
        qz = sqrt(1 - z2) * sin(phi)
        qw = sqrt(1 - z2) * cos(phi)

        @inbounds states_array[idx, 1] = 2qy * qw - 2qx * qz
        @inbounds states_array[idx, 2] = 2qy * qz + 2qx * qw
        @inbounds states_array[idx, 3] = 1 - 2qx^2 - 2qy^2
    end

    function rand_states!(states_array)
        array_shape = size(states_array)
        @assert length(array_shape) == 4
        @assert array_shape[4] == 3

        spin_size = array_shape[1:3]
        backend = KernelAbstractions.get_backend(states_array)

        kernel! = rand_states_kernel!(backend)
        kernel!(states_array, ndrange=spin_size)
    end

    @kernel function pair_energy_kernel!(
        pair_energy_array,
        center_states_array,
        interact_coeff_array,
        point_states_array,
        check_array
    )
        idx_pair = @index(Global, Cartesian)
        idx_x, idx_y, idx_t, idx_p = @index(Global, NTuple)

        if check_array[idx_x, idx_y, idx_t]
            @inbounds center_state = @view center_states_array[idx_x, idx_y, idx_t, :]
            @inbounds interact_coeff_mat = @view interact_coeff_array[idx_t, idx_p, :, :]
            @inbounds point_state = @view point_states_array[idx_pair, :]

            pair_energy = zero(eltype(pair_energy_array))
            for i = 1:3, j = 1:3
                @inbounds pair_energy += center_state[i] * interact_coeff_mat[i, j] * point_state[j]
            end

            @inbounds pair_energy_array[idx_pair] = pair_energy
        end
    end

    function pair_energy!(
        pair_energy_array,
        center_states_array,
        interact_coeff_array,
        point_states_array,
        check_array
    )
        nx_c, ny_c, nt_c, nv1_c = size(center_states_array)
        nt_i, np_i, nv1_i, nv2_i = size(interact_coeff_array)
        nx_p, ny_p, nt_p, np_p, nv2_p = size(point_states_array)
        @assert nx_c == nx_p && ny_c == ny_p && nt_c == nt_i == nt_p
        @assert np_i == np_p
        @assert nv1_c == nv1_i == nv2_i == nv2_p

        pair_size = size(pair_energy_array)
        backend = KernelAbstractions.get_backend(pair_energy_array)

        kernel! = pair_energy_kernel!(backend)
        kernel!(
            pair_energy_array,
            center_states_array,
            interact_coeff_array,
            point_states_array,
            check_array,
            ndrange=pair_size
        )
    end

    @kernel function get_point_states_kernel!(
        point_states_array,
        states_array,
        idx_array,
        check_array
    )
        idx_x, idx_y, idx_t, idx_p = @index(Global, NTuple)
        nx, ny, _, _ = size(states_array)

        if check_array[idx_x, idx_y, idx_t]
            x_diff, y_diff, point_type = @view idx_array[idx_t, idx_p, :]
            idx_x_p = mod1(idx_x + x_diff, nx)
            idx_y_p = mod1(idx_y + y_diff, ny)
            point_state = @view states_array[idx_x_p, idx_y_p, point_type, :]

            for pos_idx = 1:3
                @inbounds point_states_array[idx_x, idx_y, idx_t, idx_p, pos_idx] = point_state[pos_idx]
            end
        end
    end

    function get_point_states!(
        point_states_array,
        states_array,
        idx_array,
        check_array
    )
        pair_size = size(point_states_array)[1:4]
        backend = KernelAbstractions.get_backend(point_states_array)

        kernel! = get_point_states_kernel!(backend)
        kernel!(
            point_states_array,
            states_array,
            idx_array,
            check_array,
            ndrange=pair_size
        )
    end


    @kernel function try_flip_kernel!(
        states_array,
        rand_states_array,
        raw_energy,
        try_energy,
        check_array,
        temperature
    )
        idx_atom = @index(Global, Cartesian)

        if check_array[idx_atom]
            delta_energy = try_energy[idx_atom] - raw_energy[idx_atom]

            if delta_energy < 0 || rand() < exp(-delta_energy/temperature)
                for pos_idx = 1:3
                    @inbounds states_array[idx_atom, pos_idx] = rand_states_array[idx_atom, pos_idx]
                end
            end
        end
    end

    function try_flip!(
        states_array,
        rand_states_array,
        raw_energy,
        try_energy,
        check_array,
        temperature
    )
        atom_size = size(states_array)[1:3]
        backend = KernelAbstractions.get_backend(states_array)

        kernel! = try_flip_kernel!(backend)
        kernel!(
            states_array,
            rand_states_array,
            raw_energy,
            try_energy,
            check_array,
            temperature,
            ndrange=atom_size
        )
    end

    #TODO: MCMC
end