@testset "CrI3 Monte Carlo" begin
    mc_config = Sym4state.MC.load_config("mc_tomls/CrI3_N.TOML")

    @testset "Checkerboard" begin
        _, _, colors = Sym4state.MC.MCUtils.domain_decompose(mc_config)
        @test length(colors) == 2
    end
end

@testset "NiI2 Monte Carlo" begin
    mc_config = Sym4state.MC.load_config("mc_tomls/NiI2_NNN.TOML")

    @testset "Checkerboard" begin
        _, _, colors = Sym4state.MC.MCUtils.domain_decompose(mc_config)
        @test length(colors) == 12
    end
end