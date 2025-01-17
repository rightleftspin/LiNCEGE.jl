"""
These are some high level convenience wrappers for the NLCE process.
"""

"""
Performs simple NLCE for most use cases.

Inputs: 
      basis: basis for the unit cell, this wrapper
      assumes all atoms in the basis are the same

      primitive_vectors: primitive vectors to tile the lattice

      neighborhood: array of distances to consider as neighbors. ie. The
      first nearest neighbor will be all neighbors with distance equal 
      (up to machine precision) to the first entry of this array

      max_order: highest NLCE order to go to

Output:
    hashmap containing hashes of clusters and their corresponding multiplicities
"""
function simple_NLCE(
    basis::AbstractVector{<:AbstractVector{<:Real}},
    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    neighborhood::AbstractVector{<:Real},
    max_order::Integer,
)

    # Create the lattice
    lattice = NLCELattice(basis, primitive_vectors, neighborhood, max_order)
    # Generate clusters on the lattice
    generated_clusters = grow(lattice, max_order)
    # find all the isomorphic clusters
    iso_clusters =
        prune(isomorphic_pruning, filtering(translational_pruning, generated_clusters))

    # Account for the size of the unit cell in the pruning
    for (hash, (cluster, mult, subcluster_mult)) in iso_clusters
        if nv(cluster) > 1
            iso_clusters[hash] = (cluster, mult // length(basis), subcluster_mult)
        end
    end

    # Find all their subclusters
    subclusters = propogate(isomorphic_pruning, iso_clusters)

    # Initialize an empty output dictionary
    output_dict = Dict{AbstractNLCECluster,Vector{<:Real}}()

    # Return the final sum for all clusters
    for order = 1:max_order
        for (hash, mult) in nlce_summation(subclusters, order)
            output_dict[iso_clusters[hash][1]] =
                append!(get(output_dict, iso_clusters[hash][1], Vector{Real}()), mult)
        end
    end

    output_dict
end

"""
Performs coordinate NLCE, this is useful if your observables depend on real distance.

Inputs: 
      basis: basis for the unit cell, this wrapper
      assumes all atoms in the basis are the same

      primitive_vectors: primitive vectors to tile the lattice

      neighborhood: array of distances to consider as neighbors. ie. The
      first nearest neighbor will be all neighbors with distance equal 
      (up to machine precision) to the first entry of this array

      max_order: highest NLCE order to go to

Output:
    hashmap containing clusters and their corresponding multiplicities
"""
function coord_NLCE(
    symmetries,
    basis::AbstractVector{<:AbstractVector{<:Real}},
    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    neighborhood::AbstractVector{<:Real},
    max_order::Integer,
)

    # Create the lattice
    lattice = NLCELattice(basis, primitive_vectors, neighborhood, max_order, symmetries = symmetries)
    # Generate clusters on the lattice
    generated_clusters = grow(lattice, max_order)
    # find all the symmetrically distinct clusters
    sym_clusters =
        prune(symmetric_pruning, filtering(translational_pruning, generated_clusters))

    # Account for the size of the unit cell in the pruning
    for (hash, (cluster, mult, subcluster_mult)) in sym_clusters
        if nv(cluster) > 1
            sym_clusters[hash] = (cluster, mult // length(basis), subcluster_mult)
        end
    end

    # Find all their subclusters
    subclusters = propogate(symmetric_pruning, sym_clusters)

    # Initialize an empty output dictionary
    output_dict = Dict{AbstractNLCECluster,Vector{<:Real}}()
    cluster_hashes = Dict{AbstractNLCECluster, Integer}()

    # Return the final sum for all clusters
    for order = 1:max_order
        for (hash, mult) in nlce_summation(subclusters, order)
            output_dict[sym_clusters[hash][1]] =
                append!(get(output_dict, sym_clusters[hash][1], Vector{Real}()), mult)
        end
    end

    output_dict
end

