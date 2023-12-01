module MCExternal
    using TOML
    using Unitful: @u_str
    using UnitfulAtomic: auconvert, austrip

    using ..MCTypes
    using ..MCUtils


    const MU_B::Float32 = 0.5f0


    function array_to_vec_recursive(x)
        dims = ndims(x)
        if dims > 1
            return [
                array_to_vec_recursive(one_slice)
                for one_slice in eachslice(x, dims=dims)
            ]
        else
            return x
        end
    end

    function vec_to_array_recursion(x)
        if isbitstype(eltype(eltype(x)))
            return hcat(x...)
        else
            return stack(
                [
                    vec_to_array_recursion(item)
                    for item in x
                ]
            )
        end
    end

    # function save_config(filepath::String, )
    # TODO

    function load_config(filepath::String, T::Type=Float32)
        config = TOML.parsefile(filepath)

        # get lattice parameters
        supercell_size::Vector{Int} = config["lattice"]["size"]
        cell_mat::Array{T} = vec_to_array_recursion(config["lattice"]["cell_mat"])
        offset_mat::Array{T} = vec_to_array_recursion(config["lattice"]["offset_mat"])
        magmom_vector::Vector{T} = config["lattice"]["magmom_vector"]
        pair_mat::Array{Int} = vec_to_array_recursion(config["lattice"]["pair_mat"])
        interact_coeff_array::Array{T} = vec_to_array_recursion(config["lattice"]["interact_coeff_array"])
        interact_coeff_array = interact_coeff_array * u"meV" .|> auconvert .|> austrip

        # Dimension check
        interact_num = size(pair_mat, 2)
        @assert size(interact_coeff_array) == (3, 3, interact_num) "Mismatch between `pair_mat` and `interact_coeff_array`!"

        lattice = Lattice(
            supercell_size,
            cell_mat,
            offset_mat,
            magmom_vector,
            pair_mat,
            interact_coeff_array
        )

        # get environment parameters
        temperature_vec::Vector{T} = config["environment"]["temperature"]
        temperature_step::T = config["environment"]["temperature_step"]
        magnetic_field::Vector{T} = MU_B * config["environment"]["magnetic_field"]
        if iszero(temperature_step)
            environment_vec = [
                Environment{T}(
                    first(temperature_vec) * u"K" |> auconvert |> austrip,
                    magnetic_field
                )
            ]
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

        return lattice, environment_vec, mcmethod
    end
end