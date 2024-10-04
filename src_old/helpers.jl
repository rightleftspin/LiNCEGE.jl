"""
"""
function simple_nlce(basis, primitive_vectors, order)
    all_nlce_clusters = Dict()
    ## Generating the lattice itself
    lattice, center = generate_lattice(basis, primitive_vectors, 1, order + 3)
    for nlce_order = 1:order

        # Tile Lattice
        nlce_clusters = enumerate_subclusters(lattice, nlce_order, [center])

        # Group Clusters
        symmetric_reduction = reducing(nlce_clusters, symmetric_tagging)
        isomorphic_counts = counting(symmetric_reduction, isomorphic_tagging)

        # Tile Clusters
        subgraph_counts = Dict()
        for (key, (cluster, mult)) in isomorphic_counts
            subgraph_dict = counting(
                enumerate_subclusters(cluster, 1, Vector(1:nv(cluster))),
                isomorphic_tagging,
            )
            for subgraph_order = 2:(nv(cluster)-1)
                merge!(
                    subgraph_dict,
                    counting(
                        enumerate_subclusters(
                            cluster,
                            subgraph_order,
                            Vector(1:nv(cluster)),
                        ),
                        isomorphic_tagging,
                    ),
                )
            end
            subgraph_counts[key] = (cluster, mult, subgraph_dict)
        end
        merge!(all_nlce_clusters, subgraph_counts)
    end

    final_dict = Dict()

    for (cluster_hash, mult) in nlce_summation(all_nlce_clusters)
        final_dict[all_nlce_clusters[cluster_hash][1]] = mult
    end

    final_dict
end
