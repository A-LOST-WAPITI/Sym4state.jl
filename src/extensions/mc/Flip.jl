module MCFlip


using KernelAbstractions: @kernel, @index, @uniform, @private
using KernelAbstractions.Extras: @unroll


using ..MCTypes


export rand_states_kernel!, try_flip_kernel!


@kernel inbounds=true function rand_states_kernel!(states_array)
    idx = @index(Global, Linear)

    z2 = rand()
    phi = rand() * pi

    qx = sqrt(z2) * sin(phi)
    qy = sqrt(z2) * cos(phi)
    qz = sqrt(1 - z2) * sin(phi)
    qw = sqrt(1 - z2) * cos(phi)

    states_array[1, idx] = 2qy * qw - 2qx * qz
    states_array[2, idx] = 2qy * qz + 2qx * qw
    states_array[3, idx] = 1 - 2qx^2 - 2qy^2
end

@kernel inbounds=true function try_flip_kernel!(
    states_array,
    @Const(rand_states_array),
    @Const(pair_mat),
    @Const(interact_coeff_array),
    @Const(check_mat),
    @Const(idx_t),
    @Const(magmom_vector),
    @Const(magnetic_field),
    @Const(temperature)
)
    idx_x, idx_y = @index(Global, NTuple)
    n_x = @uniform size(states_array, 3)
    n_y = @uniform size(states_array, 4)
    n_p = @uniform size(interact_coeff_array, 3)

    # only give a try on sites those could be checked parallelly
    if check_mat[idx_x, idx_y]
        raw_state = @view states_array[:, idx_t, idx_x, idx_y]
        try_state = @view rand_states_array[:, idx_t, idx_x, idx_y]

        # init two energies
        raw_energy = zero(eltype(states_array))
        try_energy = zero(eltype(states_array))

        # energy from magnetic field
        @unroll for idx_pos = 1:3
            raw_energy -= magnetic_field[idx_pos] * raw_state[idx_pos]
            try_energy -= magnetic_field[idx_pos] * try_state[idx_pos]
        end
        raw_energy *= magmom_vector[idx_t]
        try_energy *= magmom_vector[idx_t]
        # energy from interacting
        for idx_p = 1:n_p
            pair_vec = @view pair_mat[:, idx_p]

            target_idx_x = mod1(idx_x + pair_vec[2], n_x)  # PBC
            target_idx_y = mod1(idx_y + pair_vec[3], n_y)  # PBC
            target_idx_t = pair_vec[4]

            interact_coeff_mat = @view interact_coeff_array[:, :, idx_p]
            point_state = @view states_array[:, target_idx_t, target_idx_x, target_idx_y]

            @unroll for i = 1:3
                @unroll for j = 1:3
                    raw_energy += raw_state[i] * interact_coeff_mat[i, j] * point_state[j]
                    try_energy += try_state[i] * interact_coeff_mat[i, j] * point_state[j]
                end
            end
        end

        delta_energy = try_energy - raw_energy
        if delta_energy < 0 || rand() < exp(-delta_energy/temperature)
            @unroll for idx_pos = 1:3
                raw_state[idx_pos] = try_state[idx_pos]
            end
        end
    end
end


end