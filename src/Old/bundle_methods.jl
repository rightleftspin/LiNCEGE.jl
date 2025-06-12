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

function WeakClusterExpansionBundle(
    expansion_basis::AbstractVector{<:AbstractVector{<:Real}},
    struct_per_basis::AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}},
    expansion_labels::AbstractVector{<:AbstractVector{<:Integer}},
    expansion_primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    expansion_neighbors::AbstractVector{<:Real},
    neighbors::AbstractVector{<:Real},
    max_order::Integer,
    hashing_fxn;
    basis_labels::AbstractVector{<:AbstractVector{<:Integer}} = [
        repeat([1], length(st)) for st in struct_per_basis
    ],
)
    lattice, unrotated_lattice =
        generate_primitive_lattice(expansion_primitive_vectors, max_order)

    exp_coords = add_basis_coords(expansion_basis, lattice)
    exp_sublattice_coords = add_basis_sublattice(expansion_basis, unrotated_lattice)

    unsorted_coords,
    unsorted_connections,
    unsorted_rev_connections,
    unsorted_labels,
    unsorted_translation_labels = hashing_lattice_coords(
        exp_coords,
        exp_sublattice_coords,
        struct_per_basis,
        basis_labels,
        expansion_labels,
    )

    exp_adj_list = adj_list_from_coords(exp_coords, expansion_neighbors)
    unsorted_adj_matrices = adj_matrices_weak(
        unsorted_coords,
        neighbors,
        unsorted_labels,
        unsorted_translation_labels,
        unsorted_rev_connections,
    )
    start_points = findfirst.(isapprox.(expansion_basis), (eachrow(exp_coords),))

    sort_perm_coords = sortperm(eachrow(unsorted_coords))

    coords = unsorted_coords[sort_perm_coords, :]
    labels = unsorted_labels[sort_perm_coords]
    translation_labels = unsorted_translation_labels[sort_perm_coords]
    rev_connections = unsorted_rev_connections[sort_perm_coords]
    adj_matrices = unsorted_adj_matrices[:, sort_perm_coords, sort_perm_coords]

    connections::Vector{Vector{Int}} = []
    for unsorted_connection in unsorted_connections
        temp_connection = []
        for site in unsorted_connection
            append!(temp_connection, findfirst(==(site), sort_perm_coords))
        end
        push!(connections, temp_connection)
    end

    lattice = Cluster(
        exp_adj_list,
        connections,
        adj_matrices,
        (length(unique(Iterators.flatten(basis_labels))) > 1),
        length(neighbors) != 1,
    )

    WeakClusterExpansionBundle(lattice, coords, start_points, max_order, hashing_fxn)

end


sorted_vertices = sort(unique(vcat(connections(cluster)[super_vertices]...)))
