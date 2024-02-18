module ModCore


using DocStringExtensions
using FileIO
using DataStructures: DisjointSets, IntDisjointSets, union!, num_groups, find_root!
using LinearAlgebra: norm
using ..Python
using ..Utils
using ..Types


export pre_process, post_process

"""
This function employs all the given symmetric operations in sym_op_vec to identify the analogous connections among the provided magnetic configurations in mag_struc_vec.
By doing so, it effectively reduces the total count of magnetic configurations required for calculating the interaction matrix between two atoms.

$(TYPEDSIGNATURES)
"""
function reduce_interact_mat_for_a_pair(
    mag_struc_vec::Vector{Struc},
    sym_op_vec::Vector{SymOp};
    atol=1e-2
)::IntDisjointSets{Int}
    mag_struc_count = length(mag_struc_vec)
    @assert mag_struc_count in [20, 36]

    unique_mag_struc_dict = Dict{Int, AbstractVector{Struc}}()
    fallback_ds = IntDisjointSets(mag_struc_count)
    for mag_struc in mag_struc_vec
        target_uni_num = mag_struc.uni_num
        occur_flag = false
        for (source_uni_num, one_group_struc_vec) in unique_mag_struc_dict
            for mag_struc_occur in one_group_struc_vec
                approx_flag, _ = struc_compare( # corresponding_dict is not needed here
                    mag_struc,
                    mag_struc_occur,
                    atol=atol
                )

                if approx_flag
                    union!(
                        fallback_ds,
                        source_uni_num,
                        target_uni_num
                    )
                    occur_flag = true

                    break
                end
            end

            occur_flag && break
        end

        if !occur_flag
            new_group_struc_vec = [
                sym_op * mag_struc
                for sym_op in sym_op_vec
            ]

            unique_mag_struc_dict[target_uni_num] = new_group_struc_vec
        end
    end

    return fallback_ds
end

"""
This function checks whether the size of the supercell is big enough for the given `cutoff_radius`.

$(TYPEDSIGNATURES)
"""
function supercell_check(
    py_refined_struc,  # Refined Python structure object
    mag_num_vec,  # Magnetic number vector
    mag_atom_count,  # Magnetic atom count
    cutoff_radius;  # Cutoff radius
    symprec=1e-2,  # Symmetry precision
    angle_tolerance=5.0,  # Angle tolerance
    max_supercell=10  # Maximum supercell
)
    # Local variable definitions
    local py_refined_supercell_struc, supercell_size    # Structure
    local sym_op_vec # Symmetry related
    local pair_ds, pair_relation_dict   # Pair
    local pass_flag  # Pass flag

    # Iterate over supercell ratios, from 2 to max_supercell
    for supercell_ratio = 2:max_supercell
        # Set supercell size
        supercell_size = [supercell_ratio, supercell_ratio, 1]

        # Create a supercell and copy the refined structure
        py_refined_supercell_struc = py_refined_struc.copy()
        py_refined_supercell_struc.make_supercell(supercell_size)

        # Create a symmetry analysis object for the supercell
        py_supercell_sga = py_Sga(
            py_refined_supercell_struc,
            symprec=symprec,
            angle_tolerance=angle_tolerance
        )

        # Extract a structure with only magnetism
        mag_struc = magonly(py_struc_to_struc(py_refined_supercell_struc), mag_num_vec)

        # Get the symmetry operation vector
        sym_op_vec = get_sym_op_vec(py_supercell_sga)

        # Get the center index vector and check if all equal pairs are equal
        center_idx_vec = [
            (idx - 1) * prod(supercell_size) + idx
            for idx = 1:mag_atom_count
        ]
        pass_flag, pair_ds, pair_relation_dict = equal_pair(
            mag_struc,
            center_idx_vec,
            cutoff_radius,
            sym_op_vec
        )

        # Break the loop if the condition is met
        if pass_flag
            break
        end
    end

    # Throw an error if no large enough supercell is found
    if !pass_flag
        error("Can't find large enough supercell!")
    end

    # Print information about the size of the supercell
    @info "$(supercell_size) supercell is large enough."

    # Convert Python structure to Julia structure
    struc = py_struc_to_struc(py_refined_supercell_struc)

    # Return the results
    return struc, supercell_size, sym_op_vec, pair_ds, pair_relation_dict
end


function sym4state(
    py_struc::Py,  # Python structure object
    mag_num_vec::AbstractVector{Int},  # Magnetic number vector
    cutoff_radius::T;  # Cutoff radius
    atol=1e-2,  # Absolute tolerance
    symprec=1e-2,  # Symmetry precision
    angle_tolerance=5.0,  # Angle tolerance
    max_supercell=10,  # Maximum supercell
    s_value=1.0,  # Spin value
    dump_supercell=true  # Whether to dump the supercell
) where T
    # Create a symmetry analysis object
    py_sga = py_Sga(
        py_struc,
        symprec=symprec,
        angle_tolerance=angle_tolerance
    )

    # Get the refined structure
    py_refined_struc = py_sga.get_refined_structure()

    # Get the atomic numbers
    py_num_vec = py_refined_struc.atomic_numbers
    num_vec = pyconvert(Vector{Int64}, py_num_vec)

    # Count the number of magnetic atoms
    mag_atom_count = count(∈(mag_num_vec), num_vec)

    # Get the space group number
    spg_num = pyconvert(Int, py_sga.get_space_group_number())

    # Print some information
    @info "There are $(mag_atom_count) atoms taken as magnetic in the given primitive cell."
    @info "The space group number of given structure is $(spg_num) with given `symprec`"

    # Check if the supercell is large enough
    raw_struc, supercell_size, sym_op_vec, pair_ds, pair_relation_dict = supercell_check(
        py_refined_struc,
        mag_num_vec,
        mag_atom_count,
        cutoff_radius;
        symprec=symprec,
        angle_tolerance=angle_tolerance,
        max_supercell=max_supercell
    )

    # If dump_supercell is true, dump the supercell into a file named "POSCAR"
    if dump_supercell
        isfile("POSCAR") && mv("POSCAR", "POSCAR_bak", force=true)
        py_struc.make_supercell(supercell_size)
        py_struc.to("POSCAR")
        @info "The supercell has been dumped into \"POSCAR\"."
    end

    # Return the symmetry information
    sym4state(
        raw_struc,
        supercell_size,
        mag_num_vec,
        sym_op_vec,
        pair_ds,
        pair_relation_dict,
        atol=atol,
        s_value=s_value
    )
end

function sym4state(
    raw_struc::Struc,  # Raw structure
    supercell_size::AbstractVector{Int},  # Supercell size
    mag_num_vec::AbstractVector{Int},  # Magnetic number vector
    sym_op_vec::AbstractVector{SymOp},  # Symmetry operation vector
    pair_ds::DisjointSets,  # Disjoint sets object
    pair_relation_dict::Dict;  # Dictionary
    atol::Float64=1e-2,  # Absolute tolerance
    s_value::Float64=1.0  # Spin value
)
    # Initialize the rotation symmetry flag
    rotation_sym_flag = false

    # Check if there is a rotational symmetry operation in the symmetry operation vector
    for sym_op in sym_op_vec
        if check_z_rot(sym_op)
            rotation_sym_flag = true
            break
        end
    end

    # Print a blank line
    @info ""

    # Get the number of groups in the disjoint sets object
    ngroups = num_groups(pair_ds)
    @info "There are $(ngroups) different type(s) of pairs."
    pair_vec_vec = pair_ds.revmap
    parents_vec = [
        find_root!(pair_ds, pair_vec)
        for pair_vec in pair_vec_vec
    ]
    group_parents_vec = unique(parents_vec)
    for (parent_idx, group_parent) in enumerate(group_parents_vec)
        @info "For the $(parent_idx)th group, equivalent pairs are shown as follows:"
        for (pair_vec, parent) in zip(pair_vec_vec, parents_vec)
            if parent == group_parent
                @info "$(parent) <=> $(pair_vec)"
            end
        end
    end

    # Further processing and print some information
    map_vec = Map[]
    relation_vec = CoeffMatRef[]
    map::Union{Nothing, Map} = nothing
    for (group_idx, group_parent) in enumerate(group_parents_vec)
        @info ""
        @info "For the $(group_idx)th group:"

        min_energy_num = 37 # make sure `min_energy_num` can be updated
        min_pair_vec = Int[]
        for (pair_vec, parent) in zip(pair_vec_vec, parents_vec)
            # pairs from different groups are skipped
            if parent != group_parent
                continue
            end

            temp_struc_vec = get_all_interact_struc_vec(raw_struc, mag_num_vec, pair_vec, s_value)
            mag_struc_vec = [magonly(struc, mag_num_vec) for struc in temp_struc_vec]

            temp_fallback_ds = reduce_interact_mat_for_a_pair(
                mag_struc_vec,
                sym_op_vec;
                atol=atol
            )

            energy_num = num_groups(temp_fallback_ds)
            # record the highest symmetry for now
            if energy_num < min_energy_num
                @info "Find pair $(pair_vec) with higher symmetry!"

                min_pair_vec = pair_vec
                temp_map = Map(
                    temp_fallback_ds,
                    temp_struc_vec;
                    rotation_symmetry_flag=rotation_sym_flag
                )
                min_energy_num = length(temp_map.fallback_vec)
                map = temp_map

                @info "The number of energies now is $(min_energy_num)."
            end
        end

        push!(map_vec, map)
        for (pair_vec, parent) in zip(pair_vec_vec, parents_vec)
            if parent != group_parent
                continue
            end

            pair_relation_vec = vcat(
                min_pair_vec,
                pair_vec
            )
            fixed_pair_vec = get_fixed_pair_vec(
                raw_struc,
                mag_num_vec,
                supercell_size,
                pair_vec
            )
            coeff_ref = CoeffMatRef(
                group_idx,
                fixed_pair_vec,
                pair_relation_dict[pair_relation_vec]
            )
            push!(
                relation_vec,
                coeff_ref
            )
        end
    end

    # Return the symmetry information
    return map_vec, relation_vec
end


"""
This function takes a file path, a magnetic number vector, a cutoff radius, and some optional parameters.
It performs symmetry analysis on the structure in the file and saves the results into a file named "cal.jld2".
It also prepares the inputs for VASP.

$(TYPEDSIGNATURES)
"""
function pre_process(
    filepath,  # Path of the file containing the structure
    mag_num_vec,  # Magnetic number vector
    cutoff_radius;  # Cutoff radius
    atol=1e-2,  # Absolute tolerance
    symprec=1e-2,  # Symmetry precision
    angle_tolerance=5.0,  # Angle tolerance
    max_supercell=10,  # Maximum supercell
    s_value=1.0,  # Spin value
    incar_path="./INCAR",  # Path of the INCAR file
    poscar_path="./POSCAR",  # Path of the POSCAR file
    potcar_path="./POTCAR",  # Path of the POTCAR file
    kpoints_path="./KPOINTS",  # Path of the KPOINTS file
    kwargs...  # Other optional parameters
)
    # Get the Python structure object from the file
    py_struc = get_py_struc(filepath)

    # Perform symmetry analysis
    map_vec, relation_vec = sym4state(
        py_struc,
        mag_num_vec,
        cutoff_radius;
        atol=atol,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
        max_supercell=max_supercell,
        s_value=s_value
    )

    # Print a blank line
    println()

    # Save the results into a file named "cal.jld2"
    @info "Saving the reduced map and relations into \"cal.jld2\"..."
    save(
        "cal.jld2",
        Dict(
            "map" => map_vec,
            "relation" => relation_vec
        )
    )

    # Prepare the inputs for VASP
    to_vasp_inputs(
        map_vec,
        incar_path=incar_path,
        poscar_path=poscar_path,
        potcar_path=potcar_path,
        kpoints_path=kpoints_path,
        kwargs...
    )
end

"""
This function takes a file path and some optional parameters.
It loads the symmetry information from the file and prepares the inputs for VASP.

$(TYPEDSIGNATURES)
"""
function pre_process(
    filepath;  # Path of the file containing the symmetry information
    incar_path="./INCAR",  # Path of the INCAR file
    poscar_path="./POSCAR",  # Path of the POSCAR file
    potcar_path="./POTCAR",  # Path of the POTCAR file
    kpoints_path="./KPOINTS",  # Path of the KPOINTS file
    kwargs...  # Other optional parameters
)
    # Load the symmetry information from the file
    map_vec = load(
        filepath,
        "map"
    )

    # Prepare the inputs for VASP
    to_vasp_inputs(
        map_vec,
        incar_path=incar_path,
        poscar_path=poscar_path,
        potcar_path=potcar_path,
        kpoints_path=kpoints_path,
        kwargs...
    )
end


"""
This function takes a file path.
It loads the symmetry information from the file and calculates a pair matrix and a coefficient array.

$(TYPEDSIGNATURES)
"""
function post_process(
    cal_file_path::String  # Path of the file containing the symmetry information
)
    # Get the directory of the file
    cal_dir = dirname(abspath(cal_file_path)) * "/cal/"

    # Load the symmetry information from the file
    map_vec, relation_vec = load(
        cal_file_path,
        "map",
        "relation"
    )

    # Get the pair matrix and the coefficient array
    pair_mat, coeff_array = get_pair_and_coeff(
        map_vec,
        relation_vec,
        cal_dir
    )

    # Return the pair matrix and the coefficient array
    return pair_mat, coeff_array
end



end