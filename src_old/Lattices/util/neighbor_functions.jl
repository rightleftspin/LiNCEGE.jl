function distance_neighbor_function(distances::AbstractVector{<:Real})
        function (coords::AbstractMatrix{<:Real}, sublattice_coordinates::AbstractMatrix{<:Int})
                nc = size(coords, 1)
                adj_matrix = zeros(Int, nc, nc)

                for i in 1:nc
                        for j in i+1:nc
                                dist = norm(coords[i, :] .- coords[j,:])
                                idx = findfirst(d -> isapprox(dist, d), distances)
                                if idx != nothing
                                        adj_matrix[i, j] = idx
                                        adj_matrix[j, i] = idx
                                end
                        end
                end

                adj_matrix
        end
end
