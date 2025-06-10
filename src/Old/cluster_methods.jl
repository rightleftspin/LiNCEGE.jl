begin # Complex Access functions, might take some computational effort
        function reindex_to_subcluster(
                cluster::Cluster,
                super_vertices::AbstractVector{<:Integer},
        )
                # Sort them for translational hashing
                sorted_vertices = sort(unique(vcat(connections(cluster)[super_vertices]...)))

                new_adj_list = reindex_adj_list(adj_list(cluster), super_vertices)
                new_connections =
                        reindex_connections(connections(cluster), super_vertices, sorted_vertices)
                new_adjacency_matrices = reindex_adjacency_matrices(
                        adjacency_matrices(cluster),
                        super_vertices,
                        sorted_vertices,
                )

                new_adj_list, new_connections, new_adjacency_matrices
        end

        Base.show(io::IO, cluster::Cluster) = print(
                io,
                "Cluster with $(nv(cluster)) vertices and $(length(weighted_edge_list(cluster))) bonds. Super lattice contains $(nsv(cluster)) super vertices",
        )

        # Sets default hashing of a cluster to be the translationally invariant hash
        Base.hash(cluster::Cluster, h::UInt) = hash(translational_pruning(cluster), h)
        Base.isequal(cluster1::Cluster, cluster2::Cluster) =
                (translational_pruning(cluster1) == translational_pruning(cluster2))
end
