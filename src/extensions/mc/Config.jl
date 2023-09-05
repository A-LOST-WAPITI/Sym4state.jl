module MCConfig
    using TOML
    using LinearAlgebra
    using PhysicalConstants.CODATA2018  # 参考CODATA2018标准的物理常数
    using Unitful           # 单位运算软件包
    using UnitfulAtomic     # 原子单位运算软件包


    function rotate(x::Matrix{T}, theta) where T
        r_mat::Matrix{T} = [
            cos(theta)  sin(theta)  0;
            -sin(theta) cos(theta)  0;
            0           0           1
        ]

        return r_mat * x * transpose(r_mat)
    end

    function config_cat(x::Vector{Vector{T}}) where T
        reduce(hcat, x) |> transpose |> Matrix{T}
    end

    function mev_to_au(x)
        return x * u"meV" .|> auconvert .|> austrip
    end

    function get_configs()
        μB = BohrMagneton |> auconvert |> austrip
        raw_config = TOML.parsefile("config.toml")

        configs = Dict()

        # COUPLING
        rotate_flag = raw_config["COUPLING"]["J_ROTATE"]
        all_j = raw_config["COUPLING"]["J"]
        all_a = raw_config["COUPLING"]["A"]
        all_idx = raw_config["COUPLING"]["INTERACT_INDEX"]
        temp_j = Any[]
        temp_a = Any[]
        temp_idx = Any[]
        if rotate_flag
            all_rotate_angle = raw_config["COUPLING"]["J_ROTATE_ANGLE"]
            for (rotate_angle_atom, j_atom) in zip(all_rotate_angle, all_j)
                j_mat_vec = j_atom .|> config_cat
                j_array = cat(
                    [
                        cat(
                            [rotate(j_mat, angle) for angle in deg2rad.(angle_vec)]...,
                            dims=3
                        )
                        for (j_mat, angle_vec) in zip(j_mat_vec, rotate_angle_atom)
                    ]...,
                    dims=3
                )

                push!(temp_j, j_array)
            end
        else
            for j_atom in all_j
                j_mat_vec = j_atom .|> config_cat
                j_array = cat(j_mat_vec..., dims=3)

                push!(temp_j, j_array)
            end
        end
        for a_atom in all_a
            a_mat = a_atom |> config_cat

            push!(temp_a, a_mat)
        end
        for idx_atom in all_idx
            idx_mat = vcat((idx_atom .|> config_cat)...)

            push!(temp_idx, idx_mat)
        end
        configs["j_tuple"] = tuple((temp_j .|> mev_to_au)...)
        configs["a_tuple"] = tuple((temp_a .|> mev_to_au)...)
        configs["idx_tuple"] = tuple(temp_idx...)
        # LATTICE
        configs["n_x"] = raw_config["LATTICE"]["N_LATTICE_X"]
        configs["n_y"] = raw_config["LATTICE"]["N_LATTICE_Y"]
        @assert length(temp_a) == length(temp_j) == length(temp_idx)
        configs["atom_type_num"] = length(temp_a)
        # ATOM
        m_vec = raw_config["ATOM"]["M"] * μB
        # ENV
        configs["t_max"] = raw_config["ENV"]["T_MAX"]
        configs["t_step"] = raw_config["ENV"]["T_STEP"]
        configs["t_min"] = raw_config["ENV"]["T_MIN"]
        mag_field = raw_config["ENV"]["B"]u"T" .|> auconvert .|> austrip
        configs["mag_field_tuple"] = tuple([mag_field * m for m in m_vec]...)
        # METHOD
        configs["mc_step_max"] = raw_config["METHOD"]["MC_STEP_MAX"]
        configs["measure"] = raw_config["METHOD"]["MEASURE"] |> Symbol
        configs["measure_times_max"] = raw_config["METHOD"]["MEASURE_TIMES_MAX"]
        # MISC
        configs["log"] = raw_config["MISC"]["LOG"]

        # check
        @assert length(configs["j_tuple"]) == length(configs["a_tuple"])
        return configs
    end
end