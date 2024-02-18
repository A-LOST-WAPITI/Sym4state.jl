module MCUtils
    using LazySets: convex_hull
    using DataStructures: DisjointSets, union!, find_root!
    using LinearAlgebra: norm, cross
    using ..MCTypes


    export domain_decompose


    function domain_decompose(mcconfig::MCConfig)
        pair_mat = mcconfig.pair_mat
        hull_points = convex_hull(
            eachcol(pair_mat[2:3, :]) .|> Vector |> vec
        )

        # fix the size of the lattice for being decomposible under PBC
        hull_area = Int(area_of_convex_hull(hull_points))
        size_fix_flag = !iszero(
            [
                n_lattice % hull_area
                for n_lattice in mcconfig.lattice_size
            ]
        )
        if size_fix_flag
            x_lattice = len_fix(mcconfig.lattice_size[1], hull_area)
            y_lattice = len_fix(mcconfig.lattice_size[2], hull_area)
            mcconfig = MCConfig(
                mcconfig,
                lattice_size = [x_lattice, y_lattice]
            )
            @info "Size of the lattice has been fixed to $(x_lattice)x$(y_lattice)."
        else
            x_lattice, y_lattice = mcconfig.lattice_size
        end

        equal_points_num = length(hull_points)
        equal_points_idx = [
            hull_points[idx] + hull_points[mod1(idx + 1, equal_points_num)]
            for idx = 1:equal_points_num
        ]
        color_check_mat, colors = color_check(x_lattice, y_lattice, equal_points_idx)

        return mcconfig, color_check_mat, colors
    end

    function color_check(x_lattice, y_lattice, equal_points_idx)
        idx_mat = CartesianIndices((x_lattice, y_lattice))

        check_dsu = DisjointSets{CartesianIndex{2}}(idx_mat)
        for center_idx in idx_mat
            idx_x, idx_y = center_idx.I
            for point_idx_diff in equal_points_idx
                point_idx_x = mod1(idx_x + point_idx_diff[1], x_lattice)
                point_idx_y = mod1(idx_y + point_idx_diff[2], y_lattice)
                point_idx = CartesianIndex(point_idx_x, point_idx_y)

                union!(check_dsu, center_idx, point_idx)
            end
        end

        check_mat = [
            check_dsu.intmap[find_root!(check_dsu, idx)]
            for idx in idx_mat
        ]
        colors = unique(check_mat)

        return check_mat, colors
    end


    function len_fix(n_lattice, n_window)
        round(Int, n_window * cld(n_lattice, n_window))
    end

    function area_of_convex_hull(hull_points)
        two_2_tri(x) = [x..., 0]

        start_point = hull_points[1]
        area = zero(eltype(start_point))
        for (
            point_1,
            point_2
        ) in zip(
            hull_points[2:end - 1],
            hull_points[3:end]
        )
            lin_vec_1 = two_2_tri(point_1 .- start_point)
            lin_vec_2 = two_2_tri(point_2 .- start_point)

            area += norm(cross(lin_vec_1, lin_vec_2))/2
        end

        return area
    end
end