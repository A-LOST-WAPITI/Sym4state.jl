@testset "CrI3 four-state method" begin
    py_struc = Sym4state.Python.py_Struc.from_file("CrI3/POSCAR")
    # load pymatgen and read POSCAR
    @test length(py_struc) == 8
    # get reduce resultes
    map_vec, relation_vec = Sym4state.ModCore.sym4state(
        py_struc,
        [24],
        5.0,
        dump_supercell=false
    )
    ## test `map_vec`
    @test length(map_vec) == 2
    @test length(map_vec[1].struc_vec) == 9
    @test length(map_vec[2].struc_vec) == 2
end