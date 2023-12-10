module MCCore


using ProgressMeter: Progress, next!
using KernelAbstractions: synchronize, Backend, CPU, copyto!, isgpu, partition
using KernelAbstractions: zeros as KAzeros
using Unitful: @u_str, ustrip
using UnitfulAtomic: auconvert
using Statistics: mean
using Printf: @sprintf
using LoggingExtras
using Dates


using ..MCTypes
using ..MCUtils
using ..MCFlip
using ..MCMeasure


export mcmc


function mcmc(
    mcconfig::MCConfig{T};
    backend::Backend=CPU(),
    progress_enabled::Bool=true
) where {T}
    # domain decomposition might change the lattice size
    mcconfig, color_check_mat, colors = domain_decompose(mcconfig)

    x_lattice, y_lattice = mcconfig.lattice_size
    lattice_size_tuple = (x_lattice, y_lattice)
    n_type = length(mcconfig.magmom_vector)
    n_pair = size(mcconfig.interact_coeff_array, 3)
    atom_size_tuple = (n_type, x_lattice, y_lattice)
    atom_num = prod(atom_size_tuple)
    env_num = length(mcconfig.temperature) * size(mcconfig.magnetic_field, 2)

    states_array = KAzeros(backend, T, 3, atom_size_tuple...)
    rand_states_array = KAzeros(backend, T, 3, atom_size_tuple...)
    interact_coeff_array = KAzeros(backend, T, 3, 3, n_pair)
    copyto!(backend, interact_coeff_array, mcconfig.interact_coeff_array)
    pair_mat = KAzeros(backend, Int, 4, n_pair)
    copyto!(backend, pair_mat, mcconfig.pair_mat)
    magmom_vector = KAzeros(backend, T, n_type)
    copyto!(backend, magmom_vector, mcconfig.magmom_vector)
    magnetic_field = KAzeros(backend, T, size(mcconfig.magnetic_field)...)
    copyto!(backend, magnetic_field, mcconfig.magnetic_field)
    temperature = mcconfig.temperature

    # redefine kernel functions
    if isgpu(backend)
        rand_states_kernel_with_backend! = rand_states_kernel!(backend, (32, 32), lattice_size_tuple)
        try_flip_kernel_with_backend! = try_flip_kernel!(backend, (32, 32), lattice_size_tuple)
    else
        rand_states_kernel_with_backend! = rand_states_kernel!(backend, (1024, ), lattice_size_tuple)
        try_flip_kernel_with_backend! = try_flip_kernel!(backend, (1024, ), lattice_size_tuple)
    end
    rand_states!(states_array) = rand_states_kernel_with_backend!(
        states_array,
        ndrange=lattice_size_tuple
    )
    try_flip!(x...) = try_flip_kernel_with_backend!(
        x...,
        ndrange=lattice_size_tuple
    )

    # for different atom type
    check_mat_vec = [
        begin
            checkmat_backend = KAzeros(backend, Bool, lattice_size_tuple...)
            checkmat = zeros(Bool, lattice_size_tuple...)
            checkmat[:, :] .= (color_check_mat .== color)
            copyto!(backend, checkmat_backend, checkmat)
            checkmat_backend
        end
        for color in colors
    ]
    type_idx_vec = [
        pair_mat[1, :] .== idx_t
        for idx_t = 1:n_type
    ]
    pair_mat_type_vec = [
        pair_mat[:, type_idx]
        for type_idx in type_idx_vec
    ]
    interact_coeff_array_type_vec = [
        interact_coeff_array[:, :, type_idx]
        for type_idx in type_idx_vec
    ]

    # init magnetic texture
    rand_states!(states_array)

    # loop over different environments
    norm_mean_mag_over_env = zeros(T, env_num)
    susceptibility_over_env = zeros(T, env_num)
    specific_heat_over_env = zeros(T, env_num)
    states_over_env = KAzeros(backend, T, 3, atom_size_tuple..., env_num)
    env_idx = 1
    for mag in eachcol(magnetic_field), temp in temperature
        temp_kelvin_str = @sprintf("%.4f", ustrip(auconvert(u"K", temp)))
        @info "Start equilibration progress under $(temp_kelvin_str) K."
        p = Progress(
            mcconfig.equilibration_step_num;
            showspeed=true,
            enabled=progress_enabled
        )
        for _ = 1:mcconfig.equilibration_step_num
            rand_states!(rand_states_array)
            synchronize(backend)
            for idx_t = 1:n_type
                for check_mat in check_mat_vec
                    try_flip!(
                        states_array,
                        rand_states_array,
                        pair_mat_type_vec[idx_t],
                        interact_coeff_array_type_vec[idx_t],
                        check_mat,
                        idx_t,
                        magmom_vector,
                        mag,
                        temp
                    )
                    synchronize(backend)
                end
            end

            next!(p)
        end
        @info "Start mearsuring progress under $(temp_kelvin_str) K."
        p = Progress(
            mcconfig.measuring_step_num;
            showspeed=true,
            enabled=progress_enabled
        )
        norm_mean_mag_vec = zeros(T, mcconfig.measuring_step_num)
        energy_sum_vec = zeros(T, mcconfig.measuring_step_num)
        for idx_measure = 1:mcconfig.measuring_step_num
            rand_states!(rand_states_array)
            synchronize(backend)
            for idx_t = 1:n_type
                for check_mat in check_mat_vec
                    try_flip!(
                        states_array,
                        rand_states_array,
                        pair_mat_type_vec[idx_t],
                        interact_coeff_array_type_vec[idx_t],
                        check_mat,
                        idx_t,
                        magmom_vector,
                        mag,
                        temp
                    )
                    synchronize(backend)
                end
            end

            norm_mean_mag_vec[idx_measure] = get_norm_mean_mag(states_array)
            energy_sum_vec[idx_measure] = get_energy_sum(
                states_array,
                pair_mat,
                interact_coeff_array,
                magmom_vector,
                magnetic_field
            )

            next!(p)
        end

        norm_mean_mag = mean(norm_mean_mag_vec)
        susceptibility = (mean(norm_mean_mag_vec.^2) - mean(norm_mean_mag_vec)^2)/ temp * atom_num
        specific_heat = (mean(energy_sum_vec.^2) - mean(energy_sum_vec)^2)/(atom_num * temp^2)

        norm_mean_mag_over_env[env_idx] = norm_mean_mag
        susceptibility_over_env[env_idx] = susceptibility
        specific_heat_over_env[env_idx] = specific_heat
        states_over_env[:, :, :, :, env_idx] .= states_array
        env_idx += 1
    end

    # download date from device to host
    states_over_env = Array(states_over_env)

    return states_over_env, norm_mean_mag_over_env, susceptibility_over_env, specific_heat_over_env
end


end