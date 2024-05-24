module TopoQ


using KernelAbstractions: @kernel, @index, @uniform, @localmem, @groupsize, @print, get_backend, isgpu
using KernelAbstractions.Extras: @unroll
using KernelAbstractions: zeros as KAzeros
using ..MCUtils


export q_value


# q for a single site (of the primitive cell)
function q_density(si, sj, sk)
    upper = si[1] * sj[2] * sk[3] -
        si[1] * sj[3] * sk[2] +
        si[2] * sj[3] * sk[1] -
        si[2] * sj[1] * sk[3] +
        si[3] * sj[1] * sk[2] -
        si[3] * sj[2] * sk[1]
    lower = 1
    for idx = 1:3
        lower += si[idx] * sj[idx]
        lower += si[idx] * sk[idx]
        lower += sj[idx] * sk[idx]
    end

    return 2atan(upper / lower)
end

@kernel inbounds=true function q_density_kernel!(
    density_array,
    @Const(states_array)
)
    @uniform begin
        n_x = size(states_array, 3)
        n_y = size(states_array, 4)
    end
    idx_t, idx_x, idx_y = @index(Global, NTuple)

    idx_x_p = mod1(idx_x - 1, n_x)
    idx_x_n = mod1(idx_x + 1, n_x)
    idx_y_p = mod1(idx_y - 1, n_y)
    idx_y_n = mod1(idx_y + 1, n_y)

    type_states_array = @view states_array[:, idx_t, :, :]
    @views density_array[idx_t, idx_x, idx_y] = q_density(
        type_states_array[:, idx_x, idx_y],
        type_states_array[:, idx_x_n, idx_y],
        type_states_array[:, idx_x_n, idx_y_n],
    ) + q_density(
        type_states_array[:, idx_x, idx_y],
        type_states_array[:, idx_x_p, idx_y],
        type_states_array[:, idx_x_p, idx_y_p],
    )
end

function q_value(states_array::AbstractArray{T, 4}) where T
    backend = get_backend(states_array)
    q_density_array = KAzeros(backend, T, size(states_array)[2:end])

    ndrange = size(q_density_array)
    groupsize = isgpu(backend) ? (32, ) : (1024, )
    kernel = q_density_kernel!(backend, groupsize, ndrange)
    kernel(q_density_array, states_array, ndrange=ndrange)

    # Q for different sub-lattices
    q_value_vec = sum(q_density_array, dims=(2, 3)) / (4 * pi)

    return q_value_vec
end

end