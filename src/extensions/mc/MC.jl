module MC
    include("Types.jl")
    include("Utils.jl")
    include("External.jl")
    include("Flip.jl")
    include("Measure.jl")
    include("Core.jl")


    using .MCCore
    using .MCExternal


    export mcmc, load_config
end