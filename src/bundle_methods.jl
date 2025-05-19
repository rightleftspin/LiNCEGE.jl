"""
Methods for bundles
"""

begin # Access Methods
    lattice(bundle::AbstractBundle) = bundle.lattice
    start(bundle::AbstractBundle) = bundle.start
    max_order(bundle::AbstractBundle) = bundle.max_order
    dimensions(bundle::AbstractBundle) = size(bundle.coordinates, 2)
    hashing_fxn(bundle::AbstractBundle) = bundle.hashing_fxn
    cluster_info(bundle::AbstractBundle) = bundle.cluster_info
end

begin # Writer Methods
    function set_cluster_info!(bundle::AbstractBundle, cluster_info)
        bundle.cluster_info = cluster_info
        cluster_info
    end
end

begin # Individual Cluster Methods
    # TODO: Get cluster coordinates
    function get_coordinates(bundle::AbstractBundle, cluster::Cluster)
        # act the permutation to get you to the new cluster, since
        #it may be permuted

    end
    # TODO: plot a cluster and plot a lattice here
end

begin # NLCE methods
    function initial_clusters(
        bundle::AbstractBundle,
        per_site_factor::Integer;
        single_site::Bool = false,
    )

        t_i_clusters, super_vertices = translationally_invariant_clusters(
            lattice(bundle),
            start(bundle),
            max_order(bundle),
            single_site,
            per_site_factor,
        )

        (t_i_clusters, super_vertices)
    end

    function lattice_constants!(
        bundle::AbstractBundle,
        per_site_factor::Integer,
        t_i_clusters,
        super_vertices,
    )
        cluster_info = lattice_constants(
            hashing_fxn(bundle),
            per_site_factor,
            t_i_clusters,
            super_vertices,
        )


        set_cluster_info!(bundle, cluster_info)
    end


    function subclusters!(
        bundle::AbstractBundle,
        single_site::Bool;
        hash_fxn = hashing_fxn(bundle),
    )

        for (hash, (cluster, mult, perm, svs, subclusters)) in cluster_info(bundle)
            cluster_info(bundle)[hash] = (
                cluster,
                mult,
                perm,
                svs,
                lattice_constants_only_info(
                    hash_fxn,
                    1,
                    find_subclusters(cluster, single_site)...,
                ),
            )
        end

        cluster_info(bundle)
    end

    function final_clusters(bundle::AbstractBundle, start_from_0::Bool)
        # Initialize an empty output dictionary
        output_dict = Dict{Cluster,Vector{<:Real}}()

        start = start_from_0 ? 0 : 1

        # Return the final sum for all clusters
        for order = start:max_order(bundle)
            for (hash, mult) in nlce_summation(cluster_info(bundle), order)
                output_dict[cluster_info(bundle)[hash][1]] = append!(
                    get(output_dict, cluster_info(bundle)[hash][1], Vector{Real}()),
                    mult,
                )
            end
        end

        output_dict
    end

end
