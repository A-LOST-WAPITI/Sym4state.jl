module MCMeasure
    using Statistics: mean
    using LinearAlgebra: norm
    using KernelAbstractions: @kernel, @index, get_backend
    using KernelAbstractions: zeros as KAzeros

    const MU_B::Float32 = 0.5f0

    export mag_mean, energy_mean


    function mag_mean(states_array)
        mag_mean_vec = mean(states_array, dims=(1, 2, 3))
        mag_mean = norm(mag_mean_vec)

        return mag_mean
    end

    @kernel function site_energy_kernel!(
        energy_array,
        @Const(states_array),
        @Const(point_idx_array),
        @Const(interact_coeff_array),
        @Const(magnetic_field),
    )
        idx_x, idx_y, idx_t = @index(Global, NTuple)
        n_x, n_y, _, _ = size(states_array)
        n_p = size(interact_coeff_array, 2)

        @inbounds state = @view states_array[idx_x, idx_y, idx_t, :]

        # init energy
        energy = zero(eltype(states_array))

        # energy from magnetic field
        for idx_pos = 1:3
            @inbounds energy += MU_B * magnetic_field[idx_pos] * state[idx_pos]
        end
        # energy from interacting
        for idx_p = 1:n_p
            @inbounds point_diff = @view point_idx_array[idx_t, idx_p, :]
            target_idx_x = mod1(idx_x + point_diff[1], n_x)
            target_idx_y = mod1(idx_y + point_diff[2], n_y)
            target_idx_t = point_diff[3]

            @inbounds interact_coeff_mat = @view interact_coeff_array[idx_t, idx_p, :, :]
            @inbounds point_state = @view states_array[target_idx_x, target_idx_y, target_idx_t, :]

            for i = 1:3, j = 1:3
                @inbounds energy += state[i] * interact_coeff_mat[i, j] * point_state[j]
            end
        end

        energy_array[idx_x, idx_y, idx_t] = energy
    end

    function site_energy!(
        energy_array,
        states_array,
        point_idx_array,
        interact_coeff_array,
        magnetic_field,
    )
        site_size = size(energy_array)
        backend = get_backend(energy_array)

        kernel! = site_energy_kernel!(backend)
        kernel!(
            energy_array,
            states_array,
            point_idx_array,
            interact_coeff_array,
            magnetic_field,
            ndrange=site_size
        )
    end

    function energy_mean(
        states_array::AbstractArray{T},
        point_idx_array::AbstractArray{Int},
        interact_coeff_array::AbstractArray{T},
        magnetic_field::AbstractVector{T}
    ) where T
        site_size = size(states_array)[1:3]
        backend = get_backend(states_array)
        energy_array = KAzeros(backend, T, site_size...)

        site_energy!(
            energy_array,
            states_array,
            point_idx_array,
            interact_coeff_array,
            magnetic_field
        )

        return mean(energy_array)
    end
end