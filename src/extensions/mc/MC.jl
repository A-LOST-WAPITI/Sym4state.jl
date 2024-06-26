module MC
    include("Types.jl")
    include("Utils.jl")
    include("External.jl")
    include("Flip.jl")
    include("TopoQ.jl")
    include("Core.jl")
    include("Plot.jl")


    using .MCTypes
    using .MCExternal
    using .MCCore
    using .MCPlot


    export MCConfig, load_config, save_config, mcmc, CPU, plot
end