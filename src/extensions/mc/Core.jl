module MCCore

using CUDA
using ProgressMeter: Progress, next!
using KernelAbstractions: get_backend, synchronize

using ..MCTypes
using ..MCUtils
using ..MCFlip

const use_gpu = Ref(false)
to_gpu_or_not_to_gpu(x::AbstractArray) = use_gpu[] ? CuArray(x) : x

function __init__()
    use_gpu[] = CUDA.functional()
end

function mcmc(
    lattice::Lattice{T},
    environment::Environment{T},
    mcmethod::MCMethod;
    continue_flag::Bool=false,
    prev_states_array::Union{Nothing, AbstractArray{T}}=nothing,
    progress_enabled::Bool=true
) where T
    color_check_mat, colors = domain_decompose(lattice)

    x_lattice, y_lattice = lattice.size
    n_type = length(lattice.magmom_vector)
    atom_size_tuple = (x_lattice, y_lattice, n_type)

    states_array = zeros(T, atom_size_tuple..., 3) |> to_gpu_or_not_to_gpu
    rand_states_array = zeros(T, atom_size_tuple..., 3) |> to_gpu_or_not_to_gpu
    color_check_mat = color_check_mat |> to_gpu_or_not_to_gpu
    magnetic_field = environment.magnetic_field |> to_gpu_or_not_to_gpu
    temperature = environment.temperature
    interact_coeff_array = lattice.interact_coeff_array |> to_gpu_or_not_to_gpu
    point_idx_array = lattice.point_idx_array |> to_gpu_or_not_to_gpu

    if continue_flag
        states_array .= prev_states_array
    else
        rand_states!(states_array)
    end
    check_array_mat = [
        begin
            checkarray = falses(x_lattice, y_lattice, n_type) |> to_gpu_or_not_to_gpu
            checkarray[:, :, type_idx] .= (color_check_mat .== color)
            checkarray
        end
        for color in colors, type_idx = 1:n_type
    ]
    backend = get_backend(rand_states_array)
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

    return states_array
end

end