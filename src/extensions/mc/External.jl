module MCExternal
    using TOML: parsefile

    using ..MCTypes
    using ..MCUtils


    function vec_of_vec_to_mat(
        x::AbstractVector{AbstractVector}
    )
        return permutedims( # row-major to col-major
            hcat(x...),
            (2, 1)
        )
    end

    function load_config(filepath::String)
        config = parsefile(filepath)
    end
end