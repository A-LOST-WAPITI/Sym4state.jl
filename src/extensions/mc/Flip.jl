module MCFlip


using KernelAbstractions: @kernel, @index, @uniform, @localmem, @groupsize, @print
using KernelAbstractions.Extras: @unroll


using ..MCTypes


export rand_states_kernel!, try_flip_kernel!


@inbounds function rand_state!(state::AbstractArray{T}) where T
    z2 = rand(T)
    phi = rand(T) * pi

    qx = sqrt(z2) * sin(phi)
    qy = sqrt(z2) * cos(phi)
    qz = sqrt(1 - z2) * sin(phi)
    qw = sqrt(1 - z2) * cos(phi)

    state[1] = 2qy * qw - 2qx * qz
    state[2] = 2qy * qz + 2qx * qw
    state[3] = 1 - 2qx^2 - 2qy^2
end


@kernel inbounds=true function rand_states_kernel!(states_array)
    idx = @index(Global, Cartesian)
    @uniform n_type = size(states_array, 2)

    for idx_t = 1:n_type
        rand_state!(@view states_array[:, idx_t, idx])
    end
end

@kernel inbounds=true function try_flip_kernel!(
    states_array,
    raw_energy_array,
    @Const(pair_mat),
    @Const(interact_coeff_array),
    @Const(color_idx),
    @Const(idx_t),
    @Const(magmom_vector),
    @Const(magnetic_field),
    @Const(temperature)
)
    @uniform begin
        n_x = size(states_array, 3)
        n_y = size(states_array, 4)
        n_p = size(interact_coeff_array, 3)
        group_len = prod(@groupsize())
        fixed_type = eltype(states_array)
    end

    idx_site = @index(Global, Linear)
    idx_x, idx_y = color_idx[idx_site].I
    idx_item = @index(Local, Linear)
    try_states_array = @localmem fixed_type (3, group_len)

    raw_state = @view states_array[:, idx_t, idx_x, idx_y]
    try_state = @view try_states_array[:, idx_item]

    # generate rand state
    rand_state!(try_state)

    # init two energies
    raw_energy = zero(fixed_type)
    try_energy = zero(fixed_type)

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

    raw_energy_array[idx_t, idx_x, idx_y] = raw_energy  # log site energy
    delta_energy = try_energy - raw_energy
    if delta_energy < zero(fixed_type) || rand(fixed_type) < exp(-delta_energy/temperature)
        @unroll for idx_pos = 1:3
            raw_state[idx_pos] = try_state[idx_pos]
        end
    end
end


end