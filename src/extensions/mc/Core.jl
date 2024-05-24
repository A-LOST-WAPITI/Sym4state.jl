module MCCore


using ProgressMeter: Progress, next!
using KernelAbstractions: synchronize, Backend, CPU, copyto!, isgpu, partition
using KernelAbstractions: zeros as KAzeros
using Unitful: @u_str, ustrip
using UnitfulAtomic: auconvert
using LinearAlgebra: norm
using Statistics: mean, var
using Printf: @sprintf
using LoggingExtras
using Dates


using ..MCTypes
using ..MCUtils
using ..MCFlip
using ..TopoQ


export mcmc, CPU


function mcmc(
    mcconfig::MCConfig{T};
    backend::Backend=CPU(),
    progress_enabled::Bool=true,
    log_enabled::Bool=true
) where {T}
    # domain decomposition might change the lattice size
    mcconfig, color_check_mat, colors = domain_decompose(mcconfig)

    x_lattice, y_lattice = mcconfig.lattice_size
    lattice_size_tuple = (x_lattice, y_lattice)
    color_idx_len = prod(lattice_size_tuple) / length(colors) |> Int
    n_type = length(mcconfig.magmom_vector)
    n_pair = size(mcconfig.interact_coeff_array, 3)
    atom_size_tuple = (n_type, x_lattice, y_lattice)
    atom_num = prod(atom_size_tuple)
    env_num = length(mcconfig.temperature) * size(mcconfig.magnetic_field, 2)

    # [(s_x, s_y, s_z), sites_in_primitive_cell, n_x, n_y]
    states_array = KAzeros(backend, T, 3, atom_size_tuple...)
    raw_energy_array = KAzeros(backend, T, atom_size_tuple...)
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
    groupsize = isgpu(backend) ? (32, ) : (64, )
    rand_states_kernel_with_backend! = rand_states_kernel!(backend, groupsize, lattice_size_tuple)
    try_flip_kernel_with_backend! = try_flip_kernel!(backend, groupsize, color_idx_len)
    rand_states!(states_array) = rand_states_kernel_with_backend!(
        states_array,
        ndrange=lattice_size_tuple
    )
    try_flip!(x...) = try_flip_kernel_with_backend!(
        x...,
        ndrange=color_idx_len
    )

    # for different atom type
    color_idx_vec = [
        begin
            color_idx_backend = KAzeros(backend, CartesianIndex{2}, color_idx_len)
            color_idx = findall(==(color), color_check_mat)
            copyto!(backend, color_idx_backend, color_idx)
            color_idx_backend
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
    mean_topo_q_vec_over_env = zeros(T, n_type, env_num)
    states_over_env = KAzeros(backend, T, 3, atom_size_tuple..., env_num)
    env_idx = 1
    subset_list = [
        (idx_t, color_idx)
        for idx_t = 1:n_type for color_idx in color_idx_vec
    ]
    for mag in eachcol(magnetic_field), temp in temperature
        temp_kelvin_str = @sprintf("%.4f", ustrip(auconvert(u"K", temp)))
        log_enabled && @info "Start equilibration progress under $(temp_kelvin_str) K."
        p = Progress(
            mcconfig.equilibration_step_num;
            showspeed=true,
            enabled=progress_enabled
        )
        for _ = 1:mcconfig.equilibration_step_num
            if iszero(rand((0, 1)))
                # for detailed balance https://doi.org/10.1016/j.physa.2013.08.059
                reverse!(subset_list)
            end
            for (idx_t, color_idx) in subset_list
                try_flip!(
                    states_array,
                    raw_energy_array,
                    pair_mat_type_vec[idx_t],
                    interact_coeff_array_type_vec[idx_t],
                    color_idx,
                    idx_t,
                    magmom_vector,
                    mag,
                    temp
                )
                synchronize(backend)
            end

            next!(p)
        end

        log_enabled && @info "Start mearsuring progress under $(temp_kelvin_str) K."
        p = Progress(
            mcconfig.measuring_step_num;
            showspeed=true,
            enabled=progress_enabled
        )
        norm_mean_mag_vec = zeros(T, mcconfig.measuring_step_num)
        mean_energy_vec = zeros(T, mcconfig.measuring_step_num)
        topo_q_mat = zeros(T, n_type, mcconfig.measuring_step_num)
        for idx_measure = 1:mcconfig.measuring_step_num
			for _ = 1:mcconfig.decorrelation_step_num
                if iszero(rand((0, 1)))
                    # for detailed balance https://doi.org/10.1016/j.physa.2013.08.059
                    reverse!(subset_list)
                end
                for (idx_t, color_idx) in subset_list
                    try_flip!(
                        states_array,
                        raw_energy_array,
                        pair_mat_type_vec[idx_t],
                        interact_coeff_array_type_vec[idx_t],
                        color_idx,
                        idx_t,
                        magmom_vector,
                        mag,
                        temp
                    )
                    synchronize(backend)
                end
			end

            norm_mean_mag_vec[idx_measure] = norm(mean(states_array, dims=(2, 3, 4)))
            mean_energy_vec[idx_measure] = mean(raw_energy_array) / 2
            topo_q_mat[:, idx_measure] = Array(q_value(states_array))

            next!(p)
        end

        norm_mean_mag = mean(norm_mean_mag_vec)
        mean_topo_q_vec = mean(topo_q_mat, dims=2)
        susceptibility = var(norm_mean_mag_vec; corrected=false) / temp
        specific_heat = var(mean_energy_vec; corrected=false) / temp^2
        @info @sprintf(
            "Temp = %s K: |m| = %-10.4g Chi = %-10.4g C_v = %-10.4g ",
            temp_kelvin_str, norm_mean_mag, susceptibility, specific_heat
        ) * "Q = [" * prod(x -> @sprintf("%-10.4g", x), mean_topo_q_vec) * "]"

        norm_mean_mag_over_env[env_idx] = norm_mean_mag
        susceptibility_over_env[env_idx] = susceptibility
        specific_heat_over_env[env_idx] = specific_heat
        mean_topo_q_vec_over_env[:, env_idx] .= mean_topo_q_vec
        states_over_env[:, :, :, :, env_idx] .= states_array
        env_idx += 1
    end

    # download date from device to host
    states_over_env = Array(states_over_env)

    return states_over_env, norm_mean_mag_over_env, susceptibility_over_env, specific_heat_over_env, mean_topo_q_vec_over_env
end


end
