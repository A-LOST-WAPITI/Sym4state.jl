module MCUtils
    using Printf


    function cartes_target_idx(center_idx, x_diff, y_diff, n_x, n_y; atom_type = 0)
        x, y = center_idx.I
        x_target = mod1(x + x_diff, n_x)
        y_target = mod1(y + y_diff, n_y)

        if iszero(atom_type)
            target_idx_car = CartesianIndex(x_target, y_target)
        else
            target_idx_car = CartesianIndex(x_target, y_target, atom_type)
        end

        return target_idx_car
    end

    function reducible_check(idx_atom)
        all_idx = idx_atom[:, 1:2]
        x_min, x_max = extrema(all_idx[:, 1])
        y_min, y_max = extrema(all_idx[:, 2])
        x_win = x_max - x_min + 1
        y_win = y_max - y_min + 1

        center_idx = CartesianIndex(1, 1)
        need_flag_mat = falses(x_win, y_win)
        for (x_diff, y_diff) in eachrow(all_idx)
            need_idx = cartes_target_idx(center_idx, x_diff, y_diff, x_win, y_win)
            need_flag_mat[need_idx] = true
        end

        #TODO how to properly address the check
        reducible = (x_win * y_win) - sum(need_flag_mat) >= min(x_win - 1, y_win - 1)

        if reducible
            return min(x_win, y_win), x_win, y_win
        else
            return 0, x_win, y_win
        end
    end

    function len_fix(reducible_tuple, n_x, n_y)
        para_dis_fix = lcm((reducible_tuple .|> (x -> x[1]))...)
        x_win_fix = lcm((reducible_tuple .|> (x -> x[2]))...)
        y_win_fix = lcm((reducible_tuple .|> (x -> x[3]))...)
        n_x_fix = 0
        n_y_fix = 0
        fix_flag = false

        if iszero(para_dis_fix)
            if iszero(n_x%x_win_fix) && iszero(n_y%y_win_fix)
                n_x_fix = n_x
                n_y_fix = n_y
            else
                x_time = fld(n_x, x_win_fix)
                y_time = fld(n_y, y_win_fix)

                n_x_fix = x_time * x_win_fix
                n_y_fix = y_time * y_win_fix
                fix_flag = true
            end
        else
            if iszero(n_x%para_dis_fix) && iszero(n_y%para_dis_fix)
                n_x_fix = n_x
                n_y_fix = n_y
            else
                x_time = fld(n_x, para_dis_fix)
                y_time = fld(n_y, para_dis_fix)

                n_x_fix = x_time * para_dis_fix
                n_y_fix = y_time * para_dis_fix
                fix_flag = true
            end
        end

        if fix_flag
            @warn "N_LATTICE_X has been automatically changed to $n_x_fix"
            @warn "N_LATTICE_X has been automatically changed to $n_y_fix"
        end

        return n_x_fix, n_y_fix
    end

    function dump_data(f, col_name, data)
        open(f, "w") do io
            @printf(io, "%12s\t%20s\t%20s\n", col_name...)
            for item in zip(data ...)
                @printf(io, "%12.2f\t%20.6f\t%20.6f\n", item...)
            end
        end
    end
end