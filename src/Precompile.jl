using PrecompileTools: @setup_workload, @compile_workload
using LoggingExtras
using FileIO

if Base.VERSION >= v"1.9"
    @setup_workload begin
        # Load some extra data to precompile functions
        raw_struc, supercell_size, sym_op_vec, pair_ds, pair_relation_dict = load(
            abspath(@__DIR__) * "/data/compile_input.jld2",
            "raw_struc",
            "supercell_size", 
            "sym_op_vec", 
            "pair_ds", 
            "pair_relation_dict"
        )
        log_filter = MinLevelLogger(current_logger(), Logging.Error)

        lattice, _, _ = Sym4state.MC.load_config(
            abspath(@__DIR__) * "/../test/mc_tomls/CrI3_NN.TOML"
        )
        environment = Sym4state.MC.MCTypes.Environment{Float32}(1e-5, zeros(3))
        mcmethod = Sym4state.MC.MCTypes.MCMethod(10, 1)

        with_logger(log_filter) do
            @compile_workload begin
                Sym4state.ModCore.sym4state(
                    raw_struc,
                    supercell_size,
                    [24],
                    sym_op_vec,
                    pair_ds,
                    pair_relation_dict
                )

                Sym4state.MC.mcmc(
                    lattice,
                    environment,
                    mcmethod;
                    progress_enabled=false
                )
            end
        end
    end
end