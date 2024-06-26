"""
This module defines useful functions for simplifying the four-state method.

This module exports:

$(EXPORTS)

"""
module Utils

using DocStringExtensions
using InvertedIndices: Not
using LinearAlgebra
using DelimitedFiles
using Printf
using CellListMap: neighborlist
using DataStructures: DisjointSets, find_root!, num_groups
using ..Types
using ..Python


export get_all_interact_struc_vec, to_vasp_inputs, equal_pair, get_py_struc, get_sym_op_vec, struc_compare
export linear_idx_to_vec
export py_struc_to_struc
export get_pair_and_coeff
export magonly
export check_z_rot
export get_fixed_pair_vec
export check_unit_cell


include("data/CovalentRadius.jl")

const A_COEFF_MAT::Matrix{Float64} = [
    2 1 1;
    1 2 1;
    1 1 2
]

"""
This function yields an array with size of [3, `mag_count`, 36] that comprises all 36 magnetic configurations utilized in the computation of the interaction matrix for a given number of magnetic atoms `mag_count` and the initial norm value `s_value` of [MAGMOM](https://www.vasp.at/wiki/index.php/MAGMOM). 

The resulting configurations pertain to the interaction between two atoms specified by the `target_idx_vec` vector.

$(TYPEDSIGNATURES)
"""
function mag_config(mag_count::Int, target_idx_vec::Vector{Int}, s_value::Float64)::Array{Float64, 3}
    @assert length(target_idx_vec) == 2

    axes_vec = [1, 2, 3]
    sign_vec = [1, -1]
    mag_config_vec = AbstractMatrix{Float64}[]
    idx_1, idx_2 = target_idx_vec
    axx_count = 0
    for alpha in axes_vec, beta in axes_vec
        left_axes = setdiff(axes_vec, alpha, beta)
        for sign_1 in sign_vec, sign_2 in sign_vec
            if (idx_1 == idx_2) && (alpha == beta)
                if (alpha == 1) # Axx need 4 configs
                    sign_2 = sign_1
                    # Calculing Ayy - Axx will change `left_axe` to z
                    # while calculating Azz - Axx will change `left_axe` to y
                    left_axe = iszero(axx_count % 2) ? left_axes[1] : left_axes[end]
                    axx_count += 1
                elseif iszero(sign_1 + sign_2)
                    continue
                else
                    left_axe = left_axes[end]
                end
            else    # only one axe will be left
                left_axe = left_axes[end]
            end

            mag_config = zeros(3, mag_count)
            mag_alpha = zeros(3)
            mag_beta = zeros(3)
            mag_left = zeros(3)

            mag_alpha[alpha] = sign_1
            mag_beta[beta] = sign_2
            mag_left[left_axe] = 1

            mag_config[:, idx_1] .+= mag_alpha
            mag_config[:, idx_2] .+= mag_beta
            mag_config[:, Not(target_idx_vec)] .+= mag_left

            for mag in eachcol(mag_config)
                normalize!(mag)
            end

            mag_config .*= s_value
            push!(mag_config_vec, mag_config)
        end
    end
    
    unique!(mag_config_vec)
    mag_config_array = stack(mag_config_vec)

    return mag_config_array
end

"""
This function transforms the `Structure` object from [Pymatgen](https://pymatgen.org/) into a [`Struc`](@ref) object.

The distinct identifier `uni_num` will be supplied to the `Struc` object.
The `spin_mat` variable represents the spin value assigned to each atom in the structure.

$(TYPEDSIGNATURES)
"""
function py_struc_to_struc(py_struc; uni_num::Int=1, spin_mat::Union{Nothing, AbstractMatrix}=nothing)::Struc
    py_lattice_mat = py_struc.lattice.matrix.transpose()
    py_pos_mat = py_struc.frac_coords.transpose()
    py_num_vec = py_struc.atomic_numbers

    lattice_mat = pyconvert(Matrix{Float64}, py_lattice_mat)
    pos_mat = pyconvert(Matrix{Float64}, py_pos_mat)
    # fractional coordinates should already be in [0, 1)
    # pos_mat = mod.(pos_mat, 1)
    num_vec = pyconvert(Vector{Int64}, py_num_vec)
    if isa(spin_mat, Nothing)
        spin_mat = zero(pos_mat)
    end

    struc = Struc(
        uni_num,
        lattice_mat,
        num_vec,
        pos_mat,
        spin_mat
    )

    return struc
end

"""
By utilizing the [`mag_config`](@ref),
this function generates all the [`Struc`](@ref) objects required for calculating the interaction matrix.

It takes the atomic numbers of magnetic elements stored in `mag_num_vec`,
the index of target atoms denoted by `target_idx_vec`, and the initial value `s_value`.

$(TYPEDSIGNATURES)
"""
function get_all_interact_struc_vec(
    struc::Struc,
    mag_num_vec::Vector{Int},
    target_idx_vec::Vector{Int},
    s_value::Float64
)::Vector{Struc}
    num_vec = struc.num_vec

    mag_flag_vec = [(num in mag_num_vec) for num in num_vec]
    mag_count = sum(mag_flag_vec)
    mag_config_array = mag_config(mag_count, target_idx_vec, s_value)

    struc_vec = Struc[]
    for idx in axes(mag_config_array, 3)
        spin_mat = zeros(3, length(num_vec))
        spin_mat[:, mag_flag_vec] .= mag_config_array[:, :, idx]

        struc = Struc(
            idx,
            struc.lattice_mat,
            struc.num_vec,
            struc.pos_mat,
            spin_mat
        )

        push!(struc_vec, struc)
    end

    return struc_vec
end

"""
This function examines two [`Struc`](@ref) objects and determines their approximate equivalence within a tolerance of `atol`.

It returns a `Tuple` that includes a boolean value indicating whether they are equal, 
as well as the corresponding relationships among the atoms in the respective structures.

$(TYPEDSIGNATURES)
"""
function struc_compare(x::Struc, y::Struc; atol=1e-2)::Tuple{Bool, Dict{Int, Int}}
    approx_flag = true

    corresponding_dict = Dict{Int, Int}()
    for (one_atom_idx, one_atom) in enumerate(x)
        occur_flag = false
        for (another_atom_idx, another_atom) in enumerate(y)
            if isapprox(one_atom, another_atom, atol=atol)
                occur_flag = true
                # Operation will also change the index
                # if atom i was moved to a position same as
                # raw atom j, which is to say that `y[i] = x[j]`
                corresponding_dict[another_atom_idx] = one_atom_idx
                break
            end
        end

        if !occur_flag
            approx_flag = false
            break
        end
    end

    return approx_flag, corresponding_dict
end

"""
This function examines a matrix and determines whether it represents a rotational matrix around the z-axis.

$(TYPEDSIGNATURES)
"""
function check_z_rot_mat(x::AbstractMatrix{T}; atol::T=1e-2)::Bool where T
    @assert size(x) == (3, 3)

    temp::Vector{T} = [0 ,0, 1]
    r_mat = x[1:2, 1:2]
    if !isapprox(x[:, 3], temp, atol=atol) || !isapprox(x[3, :], temp, atol=atol)
        return false
    elseif !isapprox(r_mat * transpose(r_mat), I, atol=atol) || !isapprox(det(r_mat), 1, atol=atol)
        return false
    elseif !isapprox(r_mat, I, atol=atol)
        return false
    end

    return true
end

"""
This function examines a [`SymOp`](@ref) and determines whether it contains a rotational matrix around the z-axis.

$(TYPEDSIGNATURES)
"""
check_z_rot(op::SymOp; atol=1e-2) = check_z_rot_mat(op.rot_mat, atol=atol)

"""
This function provides the pairs of atoms that are within a specified cutoff radius, 
`cutoff_radius`, and include atoms from the `center_idx_vec`.

$(TYPEDSIGNATURES)
"""
function consider_pair_vec_in_radius(
    struc::Struc,
    center_idx_vec::Vector{Int},
    cutoff_radius::Float64
)
    frac_pos_mat = struc.pos_mat
    cell_mat = struc.lattice_mat
    cart_pos_mat = cell_mat * frac_pos_mat
    neighbor_list = neighborlist(
        eachcol(cart_pos_mat),
        eachcol(cart_pos_mat),
        cutoff_radius,
        unitcell=cell_mat,
        parallel=true
    )

    consider_pair_vec_vec = Vector{Int}[]
    for pair in neighbor_list
        pair_vec = collect(pair[1:2])
        # given center idx is included
        if pair_vec[1] in center_idx_vec
            push!(
                consider_pair_vec_vec,
                pair_vec
            )
        end
    end

    return consider_pair_vec_vec
end

"""
This function computes the index of the primitive cell that encompasses the atom at index `idx`.

$(TYPEDSIGNATURES)
"""
function linear_idx_to_vec(idx::Int, supercell_size::Vector{Int}, num_pri_sites::Int)::Vector{Int}
    atom_count = prod(supercell_size) * num_pri_sites
    temp = reshape(1:atom_count, reverse(supercell_size)..., num_pri_sites)

    idx_vec = findfirst(==(idx), temp).I |> collect

    reverse!(idx_vec)

    return idx_vec
end

"""
This function calculates the corresponding pair to the `pair_vec` after a specific operation.

$(TYPEDSIGNATURES)
"""
function get_corresponding_pair_vec(
    corresponding_dict::Dict{Int, Int},
    pair_vec::AbstractVector{Int}
)::Vector{Int}
    return [
        corresponding_dict[pair_vec[1]],
        corresponding_dict[pair_vec[2]]
    ]
end

"""
This function determines all the equivalent atom pairs for a `Struc` that encloses the atom at index `center_idx_vec` within a specified cutoff radius of `cutoff_radius`. 
The symmetric operations are provided in the `sym_op_vec`.

$(TYPEDSIGNATURES)
"""
function equal_pair(
    mag_struc::Struc,
    center_idx_vec::Vector{Int},
    cutoff_radius::Float64,
    sym_op_vec::Vector{SymOp};
    atol=1e-2
)::Tuple{Bool, Union{Nothing, DisjointSets}, Union{Nothing, Dict{Vector{Int}, SymOp}}}
    # check whether unitcell is large enough
    unitcell_check_flag = check_unit_cell(mag_struc.lattice_mat, cutoff_radius)
    if !unitcell_check_flag
        return false, nothing, nothing
    end

    consider_pair_vec_vec = consider_pair_vec_in_radius(mag_struc, center_idx_vec, cutoff_radius)

    # there exists multiple pairs between center atom and one point atom
    if length(unique(consider_pair_vec_vec)) != length(consider_pair_vec_vec)
        return false, nothing, nothing
    end

    pair_ds = DisjointSets(consider_pair_vec_vec)
    pair_relation_dict = Dict{AbstractVector{Int}, SymOp}()
    for sym_op in sym_op_vec
        mag_struc_after_op = sym_op * mag_struc

        # without considering spin, those two structures should be approximate
        _, corresponding_dict = struc_compare(
            mag_struc,
            mag_struc_after_op,
            atol=atol
        )

        # find the pair after operation
        for pair_vec in consider_pair_vec_vec
            pair_vec_after_op = get_corresponding_pair_vec(
                corresponding_dict,
                pair_vec
            )
            if pair_vec_after_op in consider_pair_vec_vec
                pair_relation_vec = vcat(
                    pair_vec,
                    pair_vec_after_op
                )
                if haskey(pair_relation_dict, pair_relation_vec)
                    continue
                end

                union!(
                    pair_ds,
                    pair_vec,
                    pair_vec_after_op
                )
                pair_relation_dict[pair_relation_vec] = sym_op
            end
        end
    end

    return true, pair_ds, pair_relation_dict
end

"""
This function converts the atom pairs, 
which are identified by the indices of atoms in a supercell, 
into pairs characterized by the indices of atoms in a primitive cell along with the index of the corresponding primitive cells.

$(TYPEDSIGNATURES)
"""
function get_fixed_pair_vec(
    struc::Struc,
    mag_num_vec::Vector{Int},
    supercell_size::Vector{Int},
    pair_vec::Vector{Int}
)::Vector{Int}
    mag_struc = magonly(struc, mag_num_vec)
    mag_atom_per_cell = length(mag_struc.num_vec) ÷ prod(supercell_size)
    possibile_dis_vec = [
        (
            norm(
                mag_struc.lattice_mat * (
                    mag_struc.pos_mat[:, pair_vec[1]] - (
                        mag_struc.pos_mat[:, pair_vec[2]] .+ [diff_x, diff_y, 0]
                    )
                )
            ),
            diff_x, diff_y
        )
        for diff_x = -1:1 for diff_y in -1:1
    ]
    _, min_dis_idx = findmin(first, possibile_dis_vec)
    center_idx = linear_idx_to_vec(pair_vec[1], supercell_size, mag_atom_per_cell)
    point_idx = linear_idx_to_vec(pair_vec[2], supercell_size, mag_atom_per_cell)

    cell_idx_diff = @. (
        possibile_dis_vec[min_dis_idx][2:3] .* supercell_size[1:2] +
        point_idx[2:3] -
        center_idx[2:3]
    )
    fixed_pair_vec = vcat(center_idx[1], cell_idx_diff[1:2], point_idx[1])

    return fixed_pair_vec
end

"""
This function retrieves the [RWIGS](https://www.vasp.at/wiki/index.php/RWIGS) values for the provided path `poscar_path` of the [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR) file,
as specified in the referenced [paper](https://doi.org/10.1039/B801115J).

$(TYPEDSIGNATURES)
"""
function set_rwigs(poscar_path::String)::Vector{Float64}
    symbol_vec = split(readlines(poscar_path)[6])

    rwigs_vec = Float64[]
    for element_symbol in symbol_vec
        item = RADIUS_DICT[element_symbol]
        names = item["name"]
        radius = item["radii"]

        mult = length(names)
        if mult == 1
            radii = radius[1]
        else
            @info "There are different hybridisations of this element."
            @info "     Hybridisation: " * repr(names)
            @info "Coorespoding radii: " * repr(radius)
            @info "Type the index of the radii you want:"
            idx = parse(Int64, readline(stdin))
            
            radii = radius[idx]
        end

        @info "RWIGS of $(element_symbol) is set to $(radii)."
        push!(rwigs_vec, radii)
    end

    return rwigs_vec
end

"""
With the help of [Pymatgen](https://pymatgen.org/), this function reorganizes the `map_vec` containing [`Map`](@ref) objects and the paths of the template INCAR, POSCAR, POTCAR, and KPOINTS files into separate directories representing different magnetic configurations for subsequent DFT calculations.

Additional specific INCAR settings can be provided through the `kwargs` parameter.

$(TYPEDSIGNATURES)
"""
function to_vasp_inputs(
    map_vec::Vector{Map};
    incar_path="./INCAR",
    poscar_path="./POSCAR",
    potcar_path="./POTCAR",
    kpoints_path="./KPOINTS",
    kwargs...
)
    @assert isfile(incar_path) "Invalid `INCAR` path!"
    @assert isfile(potcar_path) "Invalid `POTCAR` path!"
    @assert isfile(kpoints_path) "Invalid `KPOINTS` path!"

    par_dir_name = "./cal/"
    if isdir(par_dir_name)
        @warn "Old input dir detected! Do you want to delete them? (Y[es]/N[o])"
        choice = readline(stdin)
        if uppercase(first(choice)) == 'Y'
            rm(par_dir_name, recursive=true, force=true)
        else
            @info "Do not delete them."
            return nothing
        end
    end

    py_incar = py_Incar.from_file(incar_path)
    rwigs_vec = set_rwigs(poscar_path)

    mkpath(par_dir_name)
    conf_dir_vec = String[]
    for (group_idx, map) in enumerate(map_vec)
        group_dir_name = par_dir_name * "group_$(group_idx)/"
        for (struc_idx, struc) in enumerate(map.struc_vec)
            conf_dir = group_dir_name * "conf_$(struc_idx)/"
            mkpath(conf_dir)

            py_magmom_list = py_np.array(transpose(struc.spin_mat)).tolist()
            update_pack_dict = Dict(
                "LNONCOLLINEAR" => true,
                "MAGMOM" => py_magmom_list,
                "M_CONSTR" => pylist(struc.spin_mat),
                "I_CONSTRAINED_M" => 1,
                "RWIGS" => pylist(rwigs_vec)
            )
            # other parameters
            for pair in kwargs
                push!(update_pack_dict, pair)
            end
            py_incar.update(update_pack_dict)

            py_incar.write_file(conf_dir * "INCAR")
            # using relative path for flexibility
            symlink(relpath(poscar_path, conf_dir), conf_dir * "POSCAR")
            symlink(relpath(potcar_path, conf_dir), conf_dir * "POTCAR")
            symlink(relpath(kpoints_path, conf_dir), conf_dir * "KPOINTS")

            push!(conf_dir_vec, conf_dir)
        end
    end

    @info "Storing path to different configuration into `cal_dir_list`. One may use SLURM's job array to calculate."
    writedlm("cal_dir_list", conf_dir_vec)
end

"""
This function using [Pymatgen](https://pymatgen.org/) to read the POSCAR file at `filepath`.

$(TYPEDSIGNATURES)
"""
get_py_struc(filepath::String) = py_Struc.from_file(filepath)

"""
This function using `pymatgen.symmetry.analyzer` from [Pymatgen](https://pymatgen.org/) to get all the symmetric operations.

$(TYPEDSIGNATURES)
"""
function get_sym_op_vec(
    py_sga
)::Vector{SymOp}
    py_sym_dict = py_sga.get_symmetry_dataset()

    lattice_mat = permutedims(
        pyconvert(Array{Float64}, py_sga._structure.lattice.matrix),
        (2, 1)
    )
    inv_lattive_mat = inv(lattice_mat)
    rot_array = permutedims(
        pyconvert(Array{Float64}, py_sym_dict["rotations"]),
        (2, 3, 1)
    )
    trans_mat = permutedims(
        pyconvert(Array{Float64}, py_sym_dict["translations"]),
        (2, 1)
    )

    op_vec = SymOp[]
    num_op = size(rot_array, 3)
    for idx = 1:num_op
        rot_mat = rot_array[:, :, idx]
        trans_vec = trans_mat[:, idx]

        for time_rev in [true, false]
            op = SymOp(
                rot_mat,
                trans_vec,
                time_rev,
                lattice_mat,
                inv_lattive_mat
            )

            push!(op_vec, op)
        end
    end

    return op_vec
end

"""
This function retrieves the energy from the OSZICAR file located at `oszicar_path`.

$(TYPEDSIGNATURES)
"""
function grep_energy(oszicar_path)
    energy = 0.0
    for line in Iterators.reverse(eachline(oszicar_path))
        m = match(r"E0", line)
        if !isnothing(m)
            energy = parse(Float64, split(line)[5])
            break
        end
    end

    # convert energy from eV to meV
    return energy * 1000
end

"""
This function computes the interaction matrix utilizing the provided [`Map`](@ref) and the energies of various configurations stored in `energy_vec`.

$(TYPEDSIGNATURES)
"""
function get_coeff_mat(map::Map, energy_vec::Dict{Int8, Float64})
    coeff_mat = zeros(3, 3)
    for i = 1:3, j = 1:3
        map_vec = map.map_mat[i, j]
        
        if map_vec[1] == 0
            coeff_mat[i, j] = 0
        else
            idx_1, idx_2, idx_3, idx_4 = map_vec
            coeff_mat[i, j] = (
                energy_vec[idx_1] -
                energy_vec[idx_2] -
                energy_vec[idx_3] +
                energy_vec[idx_4]
            )/4
        end
    end

    if map.type == 2
        coeff_mat .*= A_COEFF_MAT
    end

    return coeff_mat
end

"""
This function calculates the interaction matrix using the provided [`Map`](@ref) and the directory `cal_dir` that contains the output files of DFT calculations.

$(TYPEDSIGNATURES)
"""
function get_one_interact_coeff_mat(map::Map, cal_dir::String)
    fallback_vec = map.fallback_vec

    energy_dict = Dict{Int8, Float64}()
    for conf_idx = eachindex(fallback_vec)
        conf_dir = cal_dir * "/conf_$(conf_idx)"
        # target_conf_mag_mat = map.struc_vec[conf_idx].spin_mat
        # py_outcar = py_Outcar(conf_dir * "/OUTCAR")
        energy = grep_energy(conf_dir * "/OSZICAR")

        # target_mag_mat = mag_mat
        # target_mag_mat = round.(target_mag_mat, RoundToZero; digits=1)
        # for col in eachcol(target_mag_mat)
        #     normalize!(col)
        # end
        # if !isapprox(target_mag_mat, target_conf_mag_mat)
            # @error "MAGMOM of $(conf_dir) changed after SCF!"
        # end

        energy_dict[fallback_vec[conf_idx]] = energy
    end

    coeff_mat = get_coeff_mat(map, energy_dict)

    return coeff_mat
end

"""
This function computes all the equivalent interaction matrices within a group identified by group_id,
employing coeff_mat as the basis and considering the different relations stored in relation_vec.

$(TYPEDSIGNATURES)
"""
function get_all_interact_coeff_under_sym(
    coeff_mat::Matrix{Float64},
    group_idx::Int,
    relation_vec::Vector{CoeffMatRef}
)
    pair_vec_vec = []
    coeff_mat_vec = []
    for coeff_ref in relation_vec
        if coeff_ref.group_idx != group_idx
            continue
        end

        push!(
            pair_vec_vec,
            coeff_ref.pair_vec
        )

        push!(
            coeff_mat_vec,
            coeff_ref.op * coeff_mat
        )
    end

    pair_mat = stack(pair_vec_vec)
    coeff_array = stack(coeff_mat_vec)

    return pair_mat, coeff_array
end

"""
This function calculates all the interaction matrices and their corresponding atom pairs using the DFT calculation outputs stored in `cal_dir`.
The `map_vec` stores the relations among magnetic configurations used for calculating each interaction matrix,
while the `relation_vec` stores the relations among different interaction matrices within a group.

$(TYPEDSIGNATURES)
"""
function get_pair_and_coeff(
    map_vec::Vector{Map},
    relation_vec::Vector{CoeffMatRef},
    cal_dir::String
)::Tuple{Matrix{Int}, Array{Float64, 3}}
    group_pair_mat_vec = []
    group_coeff_array_vec = []
    for (group_idx, map) in enumerate(map_vec)
        group_dir = cal_dir * "group_$(group_idx)/"

        coeff_mat = get_one_interact_coeff_mat(map, group_dir)
        group_pair_mat, group_coeff_array = get_all_interact_coeff_under_sym(
            coeff_mat,
            group_idx,
            relation_vec
        )

        push!(
            group_pair_mat_vec,
            group_pair_mat
        )
        push!(
            group_coeff_array_vec,
            group_coeff_array
        )
    end

    pair_mat = cat(group_pair_mat_vec..., dims=2)
    coeff_array = cat(group_coeff_array_vec..., dims=3)

    return pair_mat, coeff_array
end

"""
This function filters out the non-magnetic atoms in the `struc` structure.
The atomic numbers of magnetic elements are provided in the `mag_num_vec`.

$(TYPEDSIGNATURES)
"""
function magonly(struc::Struc, mag_num_vec::Vector{Int})::Struc
    mag_flag_vec = [(num in mag_num_vec) for num in struc.num_vec]

    return Struc(
        struc.uni_num,
        struc.lattice_mat,
        struc.num_vec[mag_flag_vec],
        struc.pos_mat[:, mag_flag_vec],
        struc.spin_mat[:, mag_flag_vec]
    )
end

"""
This function checks whether the cutoff radius `cutoff` is big enough for the unit cell.
This function is designed to bypass the errors when using [CellListMap.jl](https://github.com/m3g/CellListMap.jl/) in the reference of this [file](https://github.com/m3g/CellListMap.jl/raw/3010734b6a4dd4a9366a4a1184cfffb72798b977/src/Box.jl).
$(TYPEDSIGNATURES)
"""
function check_unit_cell(unit_cell_matrix::AbstractMatrix{Float64}, cutoff::Float64)::Bool
    # from https://github.com/m3g/CellListMap.jl/raw/3010734b6a4dd4a9366a4a1184cfffb72798b977/src/Box.jl
    # remove size check and output print when check fails

    a = @view(unit_cell_matrix[:, 1])
    b = @view(unit_cell_matrix[:, 2])
    c = @view(unit_cell_matrix[:, 3])
    check = true

    bc = cross(b, c)
    bc = bc / norm(bc)
    aproj = dot(a, bc)

    ab = cross(a, b)
    ab = ab / norm(ab)
    cproj = dot(c, ab)

    ca = cross(c, a)
    ca = ca / norm(ca)
    bproj = dot(b, ca)

    if (aproj <= 2 * cutoff) || (bproj <= 2 * cutoff) || (cproj <= 2 * cutoff)
        check = false
    end

    return check
end


end
