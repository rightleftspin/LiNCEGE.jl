
function label_tiling_units(tiling_units::AbstractVector{<:AbstractMatrix{<:Real}}, primitive_vectors::AbstractMatrix{<:Real})

        pw_dir = pairwise_direction(vcat(tiling_units...))

        labels = find_labels(pw_dir, primitive_vectors)

        partition_labels(labels, tiling_units)
end

function pairwise_direction(coordinates::AbstractMatrix{<:Real})
        n_coords, dim = size(coordinates)
        pw_dir = Array{Real}(undef, n_coords, n_coords, dim)
        for i in 1:n_coords
                for j in 1:n_coords
                        pw_dir[i, j, :] = coordinates[j, :] .- coordinates[i, :]
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

        return index_matrix, unique_dirs
end

function find_labels(pw_dir::AbstractArray{<:Real,3}, primitive_vectors::AbstractMatrix{<:Real})
        n_coords, _, _ = size(pw_dir)
        labels = Vector(1:n_coords)
        prim_vector_math = transpose(primitive_vectors)

        for i in 1:n_coords
                for j in i:n_coords
                        direction = pw_dir[i, j, :]
                        int_coeffs = prim_vector_math \ direction
                        if all(isapprox.(int_coeffs, round.(int_coeffs)))
                                labels[i] = min(labels[i], labels[j])
                                labels[j] = min(labels[i], labels[j])
                        end

                end
        end

        labels
end

function partition_labels(labels::AbstractVector{<:Int}, tiling_units::AbstractVector{<:AbstractMatrix{<:Real}})
        partition_sizes = size.(tiling_units, 1)
        partitioned_labels = Vector{Vector{Int}}()
        start_idx = 1
        for partition_size in partition_sizes
                push!(partitioned_labels, labels[start_idx:partition_size+start_idx-1])
                start_idx += partition_size
        end

        partitioned_labels
end
