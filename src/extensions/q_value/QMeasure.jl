module QMeasure
    using DelimitedFiles
    using LinearAlgebra
    using HDF5
    using Printf


    # TODO: following variables as parameters
    # const all_upper_triangle_idx = [
    #     [
    #         [1, 0, 1],
    #         [1, 1, 1]
    #     ],
    # ]
    # const all_lower_triangle_idx = [
    #     [
    #         [0, -1, 1],
    #         [-1, -1, 1]
    #     ]
    # ]
    # const DATA_DIR = "data/"

    function target_idx(x, y, x_diff, y_diff, n_x, n_y)
        x_target = mod1(x + x_diff, n_x)
        y_target = mod1(y + y_diff, n_y)

        return x_target, y_target
    end

    function q_density(si, sj, sk)
        return 2atan(
            si ⋅ (sj × sk)/(
                1 + si ⋅ sj + sj ⋅ sk + sk ⋅ si
            )
        )
    end

    function julia_main()
        raw_states_array = h5read(DATA_DIR * "states.h5", "all_states")
        temp_kelvin_vec = h5read(DATA_DIR * "states.h5", "temp_vec")
        @assert size(raw_states_array, 3) == length(all_upper_triangle_idx)
        @assert size(raw_states_array, 3) == length(all_lower_triangle_idx)
        x_lim, y_lim, _ = size(raw_states_array)

        q_mat = zeros(length(temp_kelvin_vec), size(raw_states_array, 3))
        for temp_idx in eachindex(temp_kelvin_vec)
            states_array = raw_states_array[:, :, :, :, temp_idx]

            for atom_type in axes(states_array, 3)
                q_value = 0
                upper_triangle_idx = all_upper_triangle_idx[atom_type]
                lower_triangle_idx = all_lower_triangle_idx[atom_type]

                for i in axes(states_array, 1), j in axes(states_array, 2)
                    q_value += q_density(
                        states_array[i, j, atom_type, :],
                        states_array[
                            target_idx(
                                i, j,
                                upper_triangle_idx[1][1], upper_triangle_idx[1][2],
                                x_lim, y_lim
                            )..., upper_triangle_idx[1][3], 1:3
                        ],
                        states_array[
                            target_idx(
                                i, j,
                                upper_triangle_idx[2][1], upper_triangle_idx[2][2],
                                x_lim, y_lim
                            )..., upper_triangle_idx[2][3], 1:3
                        ]
                    ) + q_density(
                        states_array[i, j, atom_type, :],
                        states_array[
                            target_idx(
                                i, j,
                                lower_triangle_idx[1][1], lower_triangle_idx[1][2],
                                x_lim, y_lim
                            )..., lower_triangle_idx[1][3], 1:3
                        ],
                        states_array[
                            target_idx(
                                i, j,
                                lower_triangle_idx[2][1], lower_triangle_idx[2][2],
                                x_lim, y_lim
                            )..., lower_triangle_idx[2][3], 1:3
                        ]
                    )
                end

                q_mat[temp_idx, atom_type] = q_value/4pi
            end
        end

        writedlm(
            DATA_DIR * "q_value.dat",
            hcat(
                temp_kelvin_vec .|> (x -> @sprintf("%4.1f", x)),
                q_mat .|> (x -> @sprintf("%10.5f", x))
            )
        )
    end
end


# 如果使用脚本模式运行则自动运行`julia_main`函数
if abspath(PROGRAM_FILE) == @__FILE__
    julia_main()
end