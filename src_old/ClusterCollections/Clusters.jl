module ClusterCollections

using Base.Threads

import LINCEGE:
    _NI,
    GraphHashes.AbstractGraphHash,
    GraphHashes.IsomorphicHash,
    GraphHashes.TranslationHash,
    GraphHashes.VertexHash,
    Vertices.AbstractVertices,
    Vertices.LatticeVertices,
    Vertices.ExpansionVertices,
    Lattices.AbstractLattice,
    Lattices.AbstractSiteExpansionLattice,
    Lattices.AbstractClusterExpansionLattice,
    Lattices.max_order,
    Lattices.centers,
    Lattices.connections,
    Clusters.AbstractCluster,
    Clusters.TranslationCluster,
    Clusters.IsomorphicCluster,
    Clusters.Subgraph,
    Clusters.ghash,
    Clusters.neighbor_clusters,
    Clusters.neighbor_subgraphs,
    Clusters.merge!,
    Clusters.vertices

abstract type AbstractClusters{H<:AbstractGraphHash,C<:AbstractCluster} end

Base.length(cs::AbstractClusters) = _NI("Base.length")
Base.iterate(cs::AbstractClusters) = _NI("Base.iterate")
Base.iterate(cs::AbstractClusters, state) = _NI("Base.iterate")
Base.getindex(cs::AbstractClusters, ghash::AbstractGraphHash) = _NI("Base.getindex")
Base.setindex!(cs::AbstractClusters, cluster::AbstractCluster, ghash::AbstractGraphHash) = _NI("Base.setindex!")
Base.haskey(cs::AbstractClusters, ghash::AbstractGraphHash) = _NI("Base.haskey")

Base.in(ghash::AbstractGraphHash, cs::AbstractClusters) = haskey(cs, ghash)
Base.in(cluster::AbstractCluster, cs::AbstractClusters) = haskey(cs, ghash(cluster))
Base.show(io::IO, cs::AbstractClusters) = print(io, "Cluster Collection with $(length(cs)) clusters.")

function get_orders(cs::AbstractClusters, lattice::AbstractSiteExpansionLattice)
    result = [typeof(cs)() for i in 1:max_order(lattice)]
    for (ghash, cluster) in cs
        result[length(cluster)][ghash] = cluster
    end
    result
end

function get_orders(cs::AbstractClusters, lattice::AbstractClusterExpansionLattice)
    result = Vector{typeof(cs)}(undef, max_order(lattice))
    for (ghash, cluster) in cs
        result[length(cluster)+1][ghash] = cluster
    end
    result
end

include("TranslationClusters.jl")
include("IsomorphicClusters.jl")
include("Subgraphs.jl")

export AbstractClusters, TranslationClusters, IsomorphicClusters, Subgraphs, get_orders

end
