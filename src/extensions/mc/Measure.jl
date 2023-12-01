module MCMeasure
    using Statistics: mean
    using LinearAlgebra: norm
    using KernelAbstractions: @kernel, @index, get_backend
    using KernelAbstractions: zeros as KAzeros


    export mag_mean_norm, energy_sum


    function mag_mean_norm(states_array)
        mag_mean_vec = mean(states_array, dims=(2, 3, 4))
        result = norm(mag_mean_vec)

        return result
    end

    @kernel function site_energy_kernel!(
        energy_array,
        @Const(states_array),
        @Const(pair_mat),
        @Const(interact_coeff_array),
        @Const(magmom_vector),
        @Const(magnetic_field),
    )
        idx_site = @index(Global, Cartesian)
        idx_t, idx_x, idx_y = @index(Global, NTuple)
        _, _, n_x, n_y = size(states_array)
        n_p = size(interact_coeff_array, 3)

        @inbounds state = @view states_array[:, idx_site]

        # init energy
        energy = zero(eltype(states_array))

        # energy from magnetic field
        for idx_pos = 1:3
            @inbounds energy += magnetic_field[idx_pos] * state[idx_pos]
        end
        energy *= magmom_vector[idx_t]
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
                @inbounds energy += state[i] * interact_coeff_mat[i, j] * point_state[j]
            end
        end

        energy_array[idx_t, idx_x, idx_y] = energy
    end

    function site_energy!(
        energy_array,
        states_array,
        pair_mat,
        interact_coeff_array,
        magmom_vector,
        magnetic_field,
    )
        site_size = size(energy_array)
        backend = get_backend(energy_array)

        kernel! = site_energy_kernel!(backend)
        kernel!(
            energy_array,
            states_array,
            pair_mat,
            interact_coeff_array,
            magmom_vector,
            magnetic_field,
            ndrange=site_size
        )
    end

    function energy_sum(
        states_array::AbstractArray{T},
        pair_mat::AbstractArray{Int},
        interact_coeff_array::AbstractArray{T},
        magmom_vector::AbstractVector{T},
        magnetic_field::AbstractVector{T}
    ) where T
        site_size = size(states_array)[2:4]
        backend = get_backend(states_array)
        energy_array = KAzeros(backend, T, site_size...)

        site_energy!(
            energy_array,
            states_array,
            pair_mat,
            interact_coeff_array,
            magmom_vector,
            magnetic_field
        )

        return sum(energy_array)
    end
end