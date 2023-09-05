module MCExternal
    using TOML: parsefile
    using Rotations: RotZ
    using LinearAlgebra: transpose
    using Unitful: @u_str
    using UnitfulAtomic: auconvert, austrip

    using ..MCTypes
    using ..MCUtils


    const MU_B::Float32 = 0.5f0


    function rotate(x::AbstractMatrix{T}, theta::T) where T
        rot_mat = RotZ{T}(theta)
        return transpose(rot_mat) * x * rot_mat
    end

    function get_array_after_rotate(
        x::AbstractMatrix{T},
        theta_vec::AbstractVector{T}
    ) where T
        stack(
            [
                rotate(x, theta)
                for theta in theta_vec
            ],
            dims=1
        )
    end

    function get_interact_coeff_array_with_rotate_symmetry(
        interact_coeff_array::Array{T},
        rotate_angle
    ) where T
        temp_vec = []
        for idx_t in axes(interact_coeff_array, 1)
            temp_vec_per_type = []
            for idx_p in axes(interact_coeff_array, 2)
                interact_coeff_array_per_pair = interact_coeff_array[idx_t, idx_p, :, :]
                rotate_vec::Vector{T} = rotate_angle[idx_t][idx_p]

                push!(
                    temp_vec_per_type,
                    get_array_after_rotate(
                        interact_coeff_array_per_pair,
                        rotate_vec
                    )
                )
            end
            temp_array_per_type = cat(
                temp_vec_per_type...,
                dims=1
            )
            push!(
                temp_vec,
                temp_array_per_type
            )
        end

        return stack(
            temp_vec,
            dims=1
        )
    end

    function vec_of_vec_to_mat(
        x
    )
        return permutedims( # row-major to col-major
            hcat(x...),
            (2, 1)
        )
    end

    function vec_to_mat_recursion(x)
        if isbitstype(eltype(eltype(x)))
            return vec_of_vec_to_mat(x)
        else
            return stack(
                [
                    vec_to_mat_recursion(item)
                    for item in x
                ],
                dims=1
            )
        end
    end

    function load_config(filepath::String, T::Type=Float32)
        config = parsefile(filepath)

        # get lattice parameters
        size::Vector{Int} = config["lattice"]["size"]
        cell_mat::Array{T} = vec_of_vec_to_mat(config["lattice"]["cell_mat"])
        offset_mat::Array{T} = vec_of_vec_to_mat(config["lattice"]["offset_mat"])
        magmom_vector::Vector{T} = config["lattice"]["magmom_vector"]
        point_idx_array::Array{Int} = vec_to_mat_recursion(config["lattice"]["point_idx_array"])
        interact_coeff_array::Array{T} = vec_to_mat_recursion(config["lattice"]["interact_coeff_array"])
        rotate_symmetry::Bool = config["lattice"]["rotate_symmetry"]
        if rotate_symmetry  # if have rotation symmetry, generate array with different rotations
            rotate_angle = config["lattice"]["rotate_angle"]
            
            interact_coeff_array = get_interact_coeff_array_with_rotate_symmetry(
                interact_coeff_array,
                rotate_angle
            )
        end
        interact_coeff_array = interact_coeff_array * u"meV" .|> auconvert .|> austrip
        # TODO: Dimension check
        lattice = Lattice(
            size,
            cell_mat,
            offset_mat,
            magmom_vector,
            point_idx_array,
            interact_coeff_array
        )

        # get environment parameters
        temperature_vec::Vector{T} = config["environment"]["temperature"]
        temperature_step::T = config["environment"]["temperature_step"]
        magnetic_field::Vector{T} = MU_B * config["environment"]["magnetic_field"]
        if iszero(temperature_step)
            environment = Environment{T}(
                temperature_vec[1] * u"K" |> auconvert |> austrip,
                magnetic_field
            )
        else
            environment_vec = [
                Environment{T}(
                    temperature * u"K" |> auconvert |> austrip,
                    magnetic_field
                )
                for temperature in range(temperature_vec..., step=temperature_step)
            ]
        end

        # get mcmethod parameters
        equilibration_step_num::Int = config["mcmethod"]["equilibration_step_num"]
        measuring_step_num::Int = config["mcmethod"]["measuring_step_num"]
        mcmethod = MCMethod(
            equilibration_step_num,
            measuring_step_num
        )

        if iszero(temperature_step)
            return lattice, environment, mcmethod
        else
            return lattice, environment_vec, mcmethod
        end
    end
end