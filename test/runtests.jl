using Sym4state
using Test

@testset "Sym4state.jl" begin
    @testset "Symmetry" begin
        include("sym_tests.jl")
    end

    @testset "Monte Carlo" begin
        include("mc_tests.jl")
    end
end
