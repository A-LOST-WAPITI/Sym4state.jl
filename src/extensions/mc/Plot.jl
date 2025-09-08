module MCPlot


using LinearAlgebra
using Makie
using ColorSchemes
using Printf: @sprintf


export plot


const IDX_TYPE = Union{UnitRange{Int}, Vector{Int}}


function plot(
    states_array::AbstractArray{T, 4},  # [3, n_t, n_x, n_y]
    cell_mat::AbstractMatrix{T},    # [2, 2]
    frac_coords_pos_mat::AbstractMatrix{T}, # [2, n_t]
    target_x_idx::IDX_TYPE,
    target_y_idx::IDX_TYPE
) where T
    _, n_t, n_x, n_y = size(states_array)
    @assert size(cell_mat) == (2, 2)
    @assert size(frac_coords_pos_mat, 2) == n_t
    @assert target_x_idx ⊆ 1:n_x
    @assert target_y_idx ⊆ 1:n_y

    base_point_vec = [
        cell_mat * [i - 1, j - 1] # i for y while j for x
        for i in target_x_idx for j in target_y_idx
    ]
    cart_coords_pos_mat = cell_mat * frac_coords_pos_mat

    ps = Point3f[]
    ns = Vec3f[]
    z_vec = Float64[]
    for idx_t = 1:n_t
        ps_temp = [
            Point3f((base_point .+ cart_coords_pos_mat[:, idx_t])..., 0)
            for base_point in base_point_vec
        ]
        ns_temp = Vec3f.(
            eachslice(states_array[:, idx_t, target_x_idx, target_y_idx], dims=1)...
        )
        z_vec_temp = vec(states_array[3, idx_t, target_x_idx, target_y_idx])

        append!(ps, ps_temp)
        append!(ns, ns_temp)
        append!(z_vec, z_vec_temp)
    end

    fig = Figure(fontsize=20, backgroundcolor=:transparent)
    ax = Axis(fig[1, 1], aspect=DataAspect(), backgroundcolor=:transparent)
    hidedecorations!(ax)
    arrows3d!(
        ax,
        ps, ns,
        color=z_vec .|> (x -> get(ColorSchemes.bwr, (x + 1)/2)),
        lengthscale=1f0,
        align=:center
    )
    Colorbar(
        fig[1, 2],
        limits=(-1, 1),
        size=25,
        colormap=:bwr,
        ticklabelsize=18,
        tickformat=(X -> X .|> x -> @sprintf("%5.1f", x)),
        ticklabelalign=(:left, :center),
        height=@lift Fixed($(pixelarea(ax.scene)).widths[2])
    )
    resize_to_layout!(fig)

    return fig
end


end