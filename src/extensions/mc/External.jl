module MCExternal


using TOML
using Unitful: @u_str, ustrip
using UnitfulAtomic: auconvert, austrip
using ArgCheck
using Printf

using ..MCTypes
using ..MCUtils


export load_config, save_config


function print_array_to_vec_recursion(io::IO, x, blanks=0; no_comma_at_end=true)
    dims = ndims(x)
    prefix = "    "^blanks

    print(io, prefix * "[")
    if dims > 1
        println(io)
        for one_slice in eachslice(x, dims=dims)
            print_array_to_vec_recursion(io, one_slice, blanks + 1; no_comma_at_end=false)
        end
    else
        for item in x
            if isa(item, Integer)
                @printf(io, "%10d,", item)
            else
                @printf(io, "%10.3f,", item)
            end
        end
    end
    print(io, prefix * "]")

    if !no_comma_at_end
        println(io, ",")
    end

    return nothing
end 

function vec_to_array_recursion(x)
    if all(isbitstype, eltype.(x))
        return reduce(hcat, x)
    else
        return stack([vec_to_array_recursion(item) for item in x])
    end
end

function save_config(filepath::String, mcconfig::MCConfig)
    config_dict = Dict(
        # Lattice
        "magmom_vector" => mcconfig.magmom_vector,
        "pair_mat"      => mcconfig.pair_mat,
        "interact_coeff_array"  => ustrip.(auconvert.(u"meV", mcconfig.interact_coeff_array)),
        # Environment
        "temperature"   => ustrip.(auconvert.(u"K", mcconfig.temperature)),
        "magnetic_field"=> ustrip.(auconvert.(u"T", mcconfig.magnetic_field)),
        # Monte Carlo
        "equilibration_step_num"=> mcconfig.equilibration_step_num,
        "measuring_step_num"    => mcconfig.measuring_step_num,
        "decorrelation_step_num"=> mcconfig.decorrelation_step_num
    )

    open(filepath, "w") do io
        for key in keys(config_dict)
            print(io, key * " = ")
            if !endswith(key, "num")
                print_array_to_vec_recursion(io, config_dict[key])
            else
                @printf(io, "%10d\n", config_dict[key])
            end
            println(io)
        end
    end
end

function load_config(filepath::String, T::Type=Float32)
    config = TOML.parsefile(filepath)

    # get lattice parameters
    magmom_vector::Vector{T} = config["magmom_vector"]
    pair_mat::Array{Int} = vec_to_array_recursion(config["pair_mat"])
    interact_coeff_array::Array{T} = vec_to_array_recursion(config["interact_coeff_array"])

    # get optional lattice parameters
    lattice_size::Vector{Int} = get(config, "lattice_size", [128, 128])

    # get optional environment parameters
    temperature_step::T = get(config, "temperature_step", 0)
    temperature::Vector{T} = get(config, "temperature", zeros(T, 1)) .|> T
    if !iszero(temperature_step)
        @argcheck length(temperature) == 2 "The start and stop points should be given in `temperature` if `temperature_step` is not zero."
        temperature = temperature[1]:temperature_step:temperature[2] |> collect
    end
    magnetic_field::VecOrMat{T} = get(config, "magnetic_field", zeros(T, 3)) .|> T

    # get optional mcmethod parameters
    equilibration_step_num::Int = get(config, "equilibration_step_num", 100_000)
    measuring_step_num::Int = get(config, "measuring_step_num", 100_000)
    decorrelation_step_num::Int = get(config, "decorrelation_step_num", 10)

    mcconfig = MCConfig(
        # lattice parameters
        lattice_size=lattice_size,
        magmom_vector=magmom_vector,
        pair_mat=pair_mat,
        interact_coeff_array=interact_coeff_array,
        # environment parameters
        temperature=temperature,
        magnetic_field=magnetic_field,
        # Monte Carlo parameters
        equilibration_step_num=equilibration_step_num,
        measuring_step_num=measuring_step_num,
        decorrelation_step_num=decorrelation_step_num
    )

    return mcconfig
end


end