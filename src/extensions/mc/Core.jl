module MCCore

using ProgressMeter: Progress, next!
using KernelAbstractions: synchronize, CPU, copyto!
using KernelAbstractions: zeros as KAzeros
using Unitful: @u_str, ustrip
using UnitfulAtomic: auconvert
using Statistics: mean
using Printf: @sprintf

using ..MCTypes
using ..MCUtils
using ..MCFlip
using ..MCMeasure

function mcmc(
    lattice::Lattice{T},
    environment::Environment{T},
    mcmethod::MCMethod;
    backend=CPU(),
    continue_flag::Bool=false,
    prev_states_array::Union{Nothing, AbstractArray{T}}=nothing,
    progress_enabled::Bool=true
) where T
    color_check_mat, colors = domain_decompose(lattice)

    x_lattice, y_lattice = lattice.size
    n_type = length(lattice.magmom_vector)
    n_pair = size(lattice.interact_coeff_array, 2)
    atom_size_tuple = (x_lattice, y_lattice, n_type)

    states_array = KAzeros(backend, T, atom_size_tuple..., 3)
    rand_states_array = KAzeros(backend, T, atom_size_tuple..., 3)
    interact_coeff_array = KAzeros(backend, T, n_type, n_pair, 3, 3)
    copyto!(backend, interact_coeff_array, lattice.interact_coeff_array)
    point_idx_array = KAzeros(backend, Int, n_type, n_pair, 3)
    copyto!(backend, point_idx_array, lattice.point_idx_array)
    magnetic_field = KAzeros(backend, T, 3)
    copyto!(backend, magnetic_field, environment.magnetic_field)
    temperature = environment.temperature

    if continue_flag
        states_array .= prev_states_array
    else
        rand_states!(states_array)
    end
    check_array_mat = [
        begin
            checkarray_backend = KAzeros(backend, Bool, atom_size_tuple...)
            checkarray = zeros(Bool, atom_size_tuple...)
            checkarray[:, :, type_idx] .= (color_check_mat .== color)
            copyto!(backend, checkarray_backend, checkarray)
            checkarray_backend
        end
        for color in colors, type_idx = 1:n_type
    ]
    temp_kelvin_str = @sprintf("%.4f", ustrip(auconvert(u"K", temperature)))
    @info "Start equilibration progress under $(temp_kelvin_str) K."
    p = Progress(
        mcmethod.step_equilibration_num;
        showspeed=true,
        enabled=progress_enabled
    )
    for _ = 1:mcmethod.step_equilibration_num
        rand_states!(rand_states_array)
        synchronize(backend)
        for check_array in check_array_mat
            try_flip!(
                states_array,
                rand_states_array,
                point_idx_array,
                interact_coeff_array,
                check_array,
                magnetic_field,
                temperature
            )
            synchronize(backend)
        end

        next!(p)
    end
    @info "Start mearsuring progress under $(temp_kelvin_str) K."
    p = Progress(
        mcmethod.step_measure_num;
        showspeed=true,
        enabled=progress_enabled
    )
    mag_mean_vec = zeros(T, mcmethod.step_measure_num)
    energy_mean_vec = zeros(T, mcmethod.step_measure_num)
    for idx_measure = 1:mcmethod.step_measure_num
        rand_states!(rand_states_array)
        synchronize(backend)
        for check_array in check_array_mat
            try_flip!(
                states_array,
                rand_states_array,
                point_idx_array,
                interact_coeff_array,
                check_array,
                magnetic_field,
                temperature
            )
            synchronize(backend)
        end

        mag_mean_vec[idx_measure] = mag_mean(states_array)
        energy_mean_vec[idx_measure] = energy_mean(
            states_array,
            point_idx_array,
            interact_coeff_array,
            magnetic_field
        )

        next!(p)
    end

    mag = mean(mag_mean_vec)
    susceptibility = (mean(mag_mean_vec.^2) - mean(mag_mean_vec)^2)/temperature
    specific_heat = (mean(energy_mean_vec.^2) - mean(energy_mean_vec)^2)/temperature

    return states_array, mag, susceptibility, specific_heat
end

end