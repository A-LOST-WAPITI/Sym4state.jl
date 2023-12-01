module MCCore

using ProgressMeter: Progress, next!
using KernelAbstractions: synchronize, Backend, CPU, copyto!
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

function mcmc(
    lattice::Lattice{T},
    environment::Environment{T},
    mcmethod::MCMethod;
    backend::Backend=CPU(),
    continue_flag::Bool=false,
    prev_states_array::Union{Nothing, AbstractArray{T}}=nothing,
    progress_enabled::Bool=true
) where T
    color_check_mat, colors = domain_decompose(lattice)

    x_lattice, y_lattice = lattice.size
    n_type = length(lattice.magmom_vector)
    n_pair = size(lattice.interact_coeff_array, 3)
    atom_size_tuple = (n_type, x_lattice, y_lattice)
    atom_num = prod(atom_size_tuple)

    states_array = KAzeros(backend, T, 3, atom_size_tuple...)
    rand_states_array = KAzeros(backend, T, 3, atom_size_tuple...)
    interact_coeff_array = KAzeros(backend, T, 3, 3, n_pair)
    copyto!(backend, interact_coeff_array, lattice.interact_coeff_array)
    pair_mat = KAzeros(backend, Int, 4, n_pair)
    copyto!(backend, pair_mat, lattice.pair_mat)
    magmom_vector = KAzeros(backend, T, n_type)
    copyto!(backend, magmom_vector, lattice.magmom_vector)
    magnetic_field = KAzeros(backend, T, 3)
    copyto!(backend, magnetic_field, environment.magnetic_field)
    temperature = environment.temperature

    if continue_flag
        states_array .= prev_states_array
    else
        rand_states!(states_array)
    end
    check_array_vec = [
        begin
            checkarray_backend = KAzeros(backend, Bool, atom_size_tuple...)
            checkarray = zeros(Bool, atom_size_tuple...)
            checkarray[atom_type, :, :] .= (color_check_mat .== color)
            copyto!(backend, checkarray_backend, checkarray)
            checkarray_backend
        end
        for color in colors for atom_type = 1:n_type
    ]
    temp_kelvin_str = @sprintf("%.4f", ustrip(auconvert(u"K", temperature)))
    @info "Start equilibration progress under $(temp_kelvin_str) K."
    p = Progress(
        mcmethod.equilibration_step_num;
        showspeed=true,
        enabled=progress_enabled
    )
    for _ = 1:mcmethod.equilibration_step_num
        rand_states!(rand_states_array)
        synchronize(backend)
        for check_array in check_array_vec
            try_flip!(
                states_array,
                rand_states_array,
                pair_mat,
                interact_coeff_array,
                check_array,
                magmom_vector,
                magnetic_field,
                temperature
            )
            synchronize(backend)
        end

        next!(p)
    end
    @info "Start mearsuring progress under $(temp_kelvin_str) K."
    p = Progress(
        mcmethod.measuring_step_num;
        showspeed=true,
        enabled=progress_enabled
    )
    mag_mean_norm_vec = zeros(T, mcmethod.measuring_step_num)
    energy_sum_vec = zeros(T, mcmethod.measuring_step_num)
    for idx_measure = 1:mcmethod.measuring_step_num
        rand_states!(rand_states_array)
        synchronize(backend)
        for check_array in check_array_vec
            try_flip!(
                states_array,
                rand_states_array,
                pair_mat,
                interact_coeff_array,
                check_array,
                magmom_vector,
                magnetic_field,
                temperature
            )
            synchronize(backend)
        end

        mag_mean_norm_vec[idx_measure] = mag_mean_norm(states_array)
        energy_sum_vec[idx_measure] = energy_sum(
            states_array,
            pair_mat,
            interact_coeff_array,
            magmom_vector,
            magnetic_field
        )

        next!(p)
    end

    mag_norm = mean(mag_mean_norm_vec)
    susceptibility = (mean(mag_mean_norm_vec.^2) - mean(mag_mean_norm_vec)^2)/ temperature * atom_num
    specific_heat = (mean(energy_sum_vec.^2) - mean(energy_sum_vec)^2)/(atom_num * temperature^2)

    return states_array, mag_norm, susceptibility, specific_heat
end

function mcmc_with_environment_change(
    lattice::Lattice{T},
    environment_vec::AbstractVector{Environment{T}},
    mcmethod::MCMethod;
    backend::Backend=CPU(),
    progress_enabled::Bool=true
) where T
    states_array_vec = Vector{AbstractArray{T}}()
    mag_norm_vec = Vector{T}()
    susceptibility_vec = Vector{T}()
    specific_heat_vec = Vector{T}()

    timestamp_logger(logger) = TransformerLogger(logger) do log
        merge(log, (; message = "$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS")) $(log.message)"))
    end
    real_logger = TeeLogger(
        MinLevelLogger(timestamp_logger(FileLogger("progress.log")), Logging.Info),
        global_logger()
    )

    with_logger(real_logger) do
        states_array, mag_norm, susceptibility, specific_heat = mcmc(
            lattice,
            environment_vec[1],
            mcmethod;
            backend=backend,
            progress_enabled=progress_enabled
        )
        push!(states_array_vec, states_array)
        push!(mag_norm_vec, mag_norm)
        push!(susceptibility_vec, susceptibility)
        push!(specific_heat_vec, specific_heat)

        for environment in environment_vec[2:end]
            states_array, mag_norm, susceptibility, specific_heat = mcmc(
                lattice,
                environment,
                mcmethod;
                backend=backend,
                progress_enabled=progress_enabled,
                continue_flag=true,
                prev_states_array=states_array
            )

            push!(states_array_vec, states_array)
            push!(mag_norm_vec, mag_norm)
            push!(susceptibility_vec, susceptibility)
            push!(specific_heat_vec, specific_heat)
        end
    end

    states_change_array = stack(states_array_vec)
    return states_change_array, mag_norm_vec, susceptibility_vec, specific_heat_vec
end

end