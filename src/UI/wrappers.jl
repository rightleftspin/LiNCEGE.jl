"""
These are some high level convenience wrappers for the NLCE process.
"""

"""
Performs site expansion NLCE for most use cases.

Inputs: 
      basis: basis for the unit cell, this wrapper
      assumes all atoms in the basis are the same

      primitive_vectors: primitive vectors to tile the lattice

      neighbors: array of distances to consider as neighbors. ie. The
      first nearest neighbor will be all neighbors with distance equal 
      (up to machine precision) to the first entry of this array

      max_order: highest NLCE order to go to
      
      basis_labels: label for each site in the basis, defaults to the same for each 

Output:
    hashmap containing clusters and their corresponding multiplicities
"""
function site_expansion_NLCE(
    basis::AbstractVector{<:AbstractVector{<:Real}},
    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    neighbors::AbstractVector{<:Real},
    max_order::Integer;
    basis_labels::AbstractVector{<:Integer} = repeat([1], length(basis)),
)

    println("Setting up the lattice")
    nlce_bundle = NLCE.SiteExpansionBundle(
        basis,
        primitive_vectors,
        neighbors,
        max_order,
        NLCE.isomorphic_pruning,
        basis_labels = basis_labels,
    )

    println("Finding all translationally distinct clusters")
    t_i_clusters, super_verts =
        NLCE.initial_clusters(nlce_bundle, length(NLCE.start(nlce_bundle)))

    println("Reducing to topologically distinct clusters")
    cluster_info = NLCE.lattice_constants!(
        nlce_bundle,
        length(NLCE.start(nlce_bundle)),
        t_i_clusters,
        super_verts,
    )

    println("Finding all subclusters")
    NLCE.subclusters!(nlce_bundle, false)

    println("Returning final NLCE weights")
    (NLCE.final_clusters(nlce_bundle, false), nlce_bundle)
end

"""
Wrapper for backwards compat
"""
function simple_NLCE(
    basis::AbstractVector{<:AbstractVector{<:Real}},
    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    neighbors::AbstractVector{<:Real},
    max_order::Integer;
    basis_labels::AbstractVector{<:Integer} = repeat([1], length(basis)),
)
    site_expansion_NLCE(basis, primitive_vectors, neighbors, max_order, basis_labels)

end
