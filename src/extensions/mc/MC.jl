module MC
    include("Types.jl")
    include("Utils.jl")
    include("External.jl")
    include("Flip.jl")
    include("Core.jl")


    using .MCTypes
    using .MCExternal
    using .MCCore


    export MCConfig, load_config, save_config, mcmc, CPU
end