function get_subgraphs(c::AbstractCluster, lattice::AbstractLattice)
        if length(c) == 1
                return Set()
        end
        max_depth = length(c) - 1
        roots = [ExpansionVertices(center) for center in c.evs]
        visited = Set()

        function try_mark(cluster::ExpansionVertices)
                already = cluster in visited
                if !already
                        push!(visited, cluster)
                end

                if length(cluster) == max_depth
                        return false
                end
                !already
        end

        function dfs(cluster::ExpansionVertices)
                if !try_mark(cluster)
                        return
                end
                for ev in neighbors(lattice, cluster)
                        if ev in c.evs
                                dfs(union(cluster, ExpansionVertices(ev)))
                        end
                end
        end

        for c in roots
                dfs(c)
        end

        visited
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
