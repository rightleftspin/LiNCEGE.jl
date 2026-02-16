function pairwise_direction(coordinates::AbstractMatrix{<:Real})
    dim, n_coords = size(coordinates)
    pw_dir = Array{Real}(undef, n_coords, n_coords, dim)
    for i in 1:n_coords
        for j in 1:n_coords
            pw_dir[i, j, :] = coordinates[:, j] .- coordinates[:, i]
        end
    end

    pw_dir
end

function is_approx_in(coords::AbstractVector{<:AbstractVector{<:Real}}, coord::AbstractVector{<:Real})
    for (i, c) in enumerate(coords)
        if isapprox(c, coord, rtol=1e-8)
            return i
        end
    end

    length(coords) + 1
end

function unique_direction_indices(pw_dir::AbstractArray{<:Real,3}, bond_matrix::AbstractMatrix{Int})
    n_coords, _, _ = size(pw_dir)
    unique_dirs = Vector{Vector{Real}}()
    index_matrix = zeros(Int, n_coords, n_coords)

    for i in 1:n_coords, j in 1:n_coords
        if i == j
            index_matrix[i, j] = 0
        else
            if bond_matrix[i, j] != 0
                index_matrix[i, j] = is_approx_in(unique_dirs, pw_dir[i, j, :])
                if index_matrix[i, j] == (length(unique_dirs) + 1)
                    push!(unique_dirs, pw_dir[i, j, :])
                end
            end
        end
    end

    return index_matrix, length(unique_dirs)
end
