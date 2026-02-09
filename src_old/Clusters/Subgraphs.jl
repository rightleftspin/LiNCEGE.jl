struct Subgraph{V<:AbstractVertices,H<:VertexHash,PH<:AbstractGraphHash} <: AbstractCluster{V,H}
    vertices::V
    ghash::H
    parent_hash::PH
    parent_vertices::V
    neighbors::V
end

function Subgraph(
    vs::AbstractVertices,
    cluster::AbstractCluster{V,H},
    lattice::AbstractLattice,
) where {V,H<:IsomorphicHash}

    Subgraph(
        vs,
        VertexHash(vs),
        hashtype(cluster)(vs, lattice)[1],
        vertices(cluster),
        intersect(neighbors(lattice, vs), vertices(cluster)),
    )
end

"""
    Constructs a list of `Subgraph` objects from the neighbors of the given cluster.

    # Arguments
    - `cluster::Subgraph`: The cluster to expand around.
    - `lattice::AbstractLattice`: The lattice to expand onto.

    # Returns
    A vector of `Subgraph` objects representing the neighboring clusters.
"""
function neighbor_subgraphs(
    subgraph::Subgraph,
    lattice::AbstractLattice,
)
    ns = Subgraph[]
    for n in neighbors(subgraph)
        push!(ns, neighbor_subgraph(subgraph, eltype(subgraph)(n), lattice))
    end
    ns
end

function neighbor_subgraph(
    subgraph::Subgraph{V,H,PH},
    n::AbstractVertices,
    lattice::AbstractLattice,
) where {V,H,PH<:IsomorphicHash}
    new_vs = union(vertices(subgraph), n)

    Subgraph(
        new_vs,
        VertexHash(new_vs),
        parent_hash_type(subgraph)(new_vs, lattice)[1],
        parent_vertices(subgraph),
        setdiff(union(neighbors(subgraph), intersect(neighbors(lattice, n), parent_vertices(subgraph))), new_vs),
    )
end

vertices(cluster::Subgraph) = cluster.vertices
ghash(cluster::Subgraph) = cluster.ghash
lattice_constant(cluster::Subgraph) = 1

parent_hash(cluster::Subgraph) = cluster.parent_hash
parent_hash_type(cluster::Subgraph{V,H,PH}) where {V,H,PH} = PH
parent_vertices(cluster::Subgraph) = cluster.parent_vertices
neighbors(cluster::Subgraph) = cluster.neighbors
