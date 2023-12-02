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
            end
        end
    end
end