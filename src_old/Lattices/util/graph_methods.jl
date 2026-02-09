function adj_mat_to_adj_list(adj_matrix::AbstractMatrix{<:Int})
    n_coords = size(adj_matrix, 1)
    adj_list = Vector{LatticeVertices}(undef, n_coords)
    for i in 1:n_coords
        neighbors = LatticeVertices()
        for j in 1:n_coords
            if adj_matrix[i, j] > 0
                neighbors = union(neighbors, LatticeVertices(j))
            end
        end
        adj_list[i] = neighbors
    end
    adj_list
end

function adj_mat_to_adj_list_exp(adj_matrix::AbstractMatrix{<:Int})
    n_coords = size(adj_matrix, 1)
    adj_list = Vector{ExpansionVertices}(undef, n_coords)
    for i in 1:n_coords
        neighbors = ExpansionVertices()
        for j in 1:n_coords
            if adj_matrix[i, j] > 0
                neighbors = union(neighbors, ExpansionVertices(j))
            end
        end
        adj_list[i] = neighbors
    end
    adj_list
end

function pairwise_distance_mat_to_adj_mat(pw_dist::AbstractMatrix{<:Real}, distances::AbstractVector{<:Real})
    n_coords, _ = size(pw_dist)
    index_matrix = Matrix{Int}(undef, size(pw_dist))
    dist_map = Dict(d => i for (i, d) in enumerate(distances))
    for i in 1:n_coords
        for j in 1:n_coords
            index_matrix[i, j] = get(dist_map, pw_dist[i, j], 0)
        end
    end

    index_matrix
end

function adj_mat_to_edge_list(adj_matrix::AbstractMatrix{<:Real})
    edge_list = []

    n = size(adj_matrix, 1)
    for i in 1:n
        for j in i+1:n
            weight = adj_matrix[i, j]
            if weight != 0
                push!(edge_list, (i, j, weight))
            end
        end
    end

    edge_list
end
