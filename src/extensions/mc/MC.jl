module MC
    include("Types.jl")
    include("Utils.jl")
    include("External.jl")
    include("Flip.jl")
    include("Measure.jl")
    include("Core.jl")

    using .MCCore: mcmc, mcmc_with_environment_change
    using .MCExternal: load_config

    export mcmc, mcmc_with_environment_change, load_config
end