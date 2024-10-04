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
function simple_NLCE(basis::AbstractVector{<:AbstractVector{T}}, primitive_vectors::AbstractVector{<:AbstractVector{T}}, neighborhood::AbstractVector{T}, max_order::Integer) where T<:Real

    # Create the lattice
    lattice = NLCELattice(basis, primitive_vec, neighborhood, max_order)
    # Generate clusters on the lattice
    generated_clusters = grow(lattice, max_order)
    # find all the isomorphic clusters
    iso_clusters = prune(isomorphic_pruning, filtering(symmetric_pruning, generated_clusters))
    # Find all their subclusters
    subclusters = propogate(isomorphic_pruning, iso_clusters) 

    # Initialize an empty output dictionary
    output_dict = Dict{AbstractNLCECluster, Vector{<:Integer}}()

    # Return the final sum for all clusters
    for order in 1:max_order
        for (hash, mult) in nlce_summation(subclusters, order)
            output_dict[iso_clusters[hash][1]] = append!(get(output_dict, iso_clusters[hash][1], Vector{Int}()), mult)
        end
    end

    output_dict

end



