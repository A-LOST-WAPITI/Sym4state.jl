module MCExternal


using TOML
using Unitful: @u_str
using UnitfulAtomic: auconvert, austrip
using ArgCheck

using ..MCTypes
using ..MCUtils


export load_config


const MU_B::Rational{Int} = 1//2


### This function converts a multi-dimensional array into a list of its slices.
### It calls itself recursively to handle arrays of any dimension.
### 
### @param x Multi-dimensional array to be converted.
### @return If the array is of dimension greater than 1, it 
### splits it on the last dimension,
### convert each slice recursively to a list, and finally return the list of 
### all slices. If the array dimension is 1, it simply returns the array.
###
function array_to_vec_recursive(x)
    dims = ndims(x) # Get the number of dimensions of the array x

    # If the array has more than one dimension
    if dims > 1
        # Recursively process each slice of the array by passing each slice to
        # the same function. The `eachslice()` function returns an iterable 
        # over the array slices in the specified dimension.
        return [
            array_to_vec_recursive(one_slice)
            for one_slice in eachslice(x, dims=dims)
        ]
    else
        # If the array has only one dimension, return the array itself
        return x
    end
end 


function vec_to_array_recursion(x)
    if all(isbitstype, eltype.(x))
        return reduce(hcat, x)
    else
        return stack([vec_to_array_recursion(item) for item in x])
    end
end


function load_config(filepath::String, T::Type=Float32)
    config = TOML.parsefile(filepath)

    # get lattice parameters
    magmom_vector::Vector{T} = config["magmom_vector"]
    pair_mat::Array{Int} = vec_to_array_recursion(config["pair_mat"])
    interact_coeff_array::Array{T} = vec_to_array_recursion(config["interact_coeff_array"])
    interact_coeff_array = @. interact_coeff_array * u"meV" |> auconvert |> austrip

    # using parameters to construct a `MCConfig`
    mcconfig = MCConfig(
        magmom_vector,
        pair_mat,
        interact_coeff_array
    )

    # get optional lattice parameters
    lattice_size::Vector{Int} = get(config, "lattice_size", [128, 128])

    # get optional environment parameters
    temperature_step::T = get(config, "temperature_step", 0)
    temperature::Vector{T} = get(config, "temperature", zeros(T, 1))
    if !iszero(temperature_step)
        @argcheck length(temperature) == 2 "The start and stop points should be given in `temperature` if `temperature_step` is not zero."
        temperature = temperature[1]:temperature_step:temperature[2] |> collect
    end
    temperature = @. temperature * u"K" |> auconvert |> austrip
    magnetic_field::VecOrMat{T} = MU_B * get(config, "magnetic_field", zeros(T, 3)) .|> T
    magnetic_field = @. magnetic_field * u"T" |> auconvert |> austrip |> T

    # get optional mcmethod parameters
    equilibration_step_num::Int = get(config, "equilibration_step_num", 100_000)
    measuring_step_num::Int = get(config, "measuring_step_num", 100_000)

    mcconfig = MCConfig(
        mcconfig,
        # update lattice parameters
        lattice_size=lattice_size,
        # update environment parameters
        temperature=temperature,
        magnetic_field=magnetic_field,
        # update Monte Carlo parameters
        equilibration_step_num=equilibration_step_num,
        measuring_step_num=measuring_step_num
    )

    return mcconfig
end


end