"""
    TranslationCluster{V<:AbstractVertices} <: AbstractCluster

A data structure that represents a cluster on a given lattice.

# Fields
- `vertices::V`: The vertices of the cluster.
- `neighbors::V`: The vertices neighboring the cluster.
"""
struct TranslationCluster{V<:AbstractVertices,H<:TranslationHash} <: AbstractCluster{V,H}
    vertices::V
    ghash::H
    neighbors::V
end

function TranslationCluster(
    vs::AbstractVertices,
    lattice::AbstractLattice,
)

    TranslationCluster(
        vs,
        TranslationHash(vs, lattice),
        neighbors(lattice, vs),
    )
end

"""
    Constructs a list of `TranslationCluster` objects from the neighbors of the given cluster.

    # Arguments
    - `cluster::TranslationCluster`: The cluster to expand around.
    - `lattice::AbstractLattice`: The lattice to expand onto.

    # Returns
    A vector of `TranslationCluster` objects representing the neighboring clusters.
"""
function neighbor_clusters(
    cluster::TranslationCluster,
    lattice::AbstractLattice,
)
    ns = TranslationCluster[]
    for n in neighbors(cluster)
        push!(ns, neighbor_cluster(cluster, eltype(cluster)(n), lattice))
    end
    ns
end

function neighbor_cluster(
    cluster::TranslationCluster,
    n::AbstractVertices,
    lattice::AbstractLattice,
)
    new_vs = union(vertices(cluster), n)

    TranslationCluster(
        new_vs,
        TranslationHash(new_vs, lattice),
        setdiff(union(neighbors(cluster), neighbors(lattice, n)), new_vs),
    )
end

vertices(cluster::TranslationCluster) = cluster.vertices
ghash(cluster::TranslationCluster) = cluster.ghash
lattice_constant(cluster::TranslationCluster) = 1

neighbors(cluster::TranslationCluster) = cluster.neighbors
