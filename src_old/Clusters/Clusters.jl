module Clusters

import LINCEGE:
    Vertices.AbstractVertices,
    Lattices.AbstractLattice,
    Lattices.neighbors,
    GraphHashes.AbstractGraphHash,
    GraphHashes.TranslationHash,
    GraphHashes.IsomorphicHash,
    GraphHashes.VertexHash,
    GraphHashes.AbstractPermutation,
    GraphHashes.EmptyPermutation,
    GraphHashes.IsomorphicPermutation,
    _NI

abstract type AbstractCluster{V<:AbstractVertices,H<:AbstractGraphHash} end

vertices(cluster::AbstractCluster) = _NI("vertices")
ghash(cluster::AbstractCluster) = _NI("ghash")
lattice_constant(cluster::AbstractCluster) = _NI("lattice_constant")

Base.eltype(cluster::AbstractCluster{V,H}) where {V,H} = V
Base.length(cluster::AbstractCluster) = length(vertices(cluster))
Base.show(io::IO, cluster::AbstractCluster) = print(io, "Cluster with $(vertices(cluster)) and Lattice Constant $(lattice_constant(cluster))")
hashtype(cluster::AbstractCluster{V,H}) where {V,H} = H

include("TranslationClusters.jl")
include("IsomorphicClusters.jl")
include("Subgraphs.jl")

export AbstractCluster,
    TranslationCluster,
    IsomorphicCluster,
    Subgraph,
    vertices,
    ghash,
    lattice_constant
end
