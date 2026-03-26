module Clusters

using Base.Threads
using LinearAlgebra
using NautyGraphs

import LINCEGE:
        _NI,
        Vertices.ExpansionVertices,
        Lattices.centers,
        Lattices.max_order,
        Lattices.neighbors,
        Lattices.n_unique_sites,
        Lattices.get_single_site_subgraphs,
        Lattices.AbstractLattice,
        Lattices.AbstractInfiniteLattice,
        Hashers.AbstractHasher,
        Hashers.TranslationHasher,
        Hashers.IsomorphicHasher,
        Hashers.ghash

abstract type AbstractCluster end
abstract type AbstractClusterSet{C<:AbstractCluster,H<:AbstractHasher} end

# Cluster Methods
Base.length(c::AbstractCluster) = _NI("Base.length")
Base.hash(c::AbstractCluster, h::UInt) = _NI("Base.hash")

Base.isequal(c1::C, c2::C) where {C<:AbstractCluster} = c1 == c2
Base.:(==)(c1::C, c2::C) where {C<:AbstractCluster} = (hash(c1) == hash(c2))

# Cluster Set Methods
Base.length(cs::AbstractClusterSet)::Int = _NI("Base.length")
Base.in(cluster::C, cs::AbstractClusterSet{C,H}) where {C<:AbstractCluster,H} = _NI("Base.in")
Base.iterate(cs::AbstractClusterSet) = _NI("Base.iterate")
Base.iterate(cs::AbstractClusterSet, state) = _NI("Base.iterate")
Base.push!(cs::AbstractClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = _NI("Base.push!")
Base.pop!(cs::AbstractClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = _NI("Base.pop!")
ghash(cs::AbstractClusterSet, c::AbstractCluster) = _NI("ghash")
ghash(cs::AbstractClusterSet, evs::ExpansionVertices) = _NI("ghash")

include("util.jl")
include("Cluster.jl")
include("ClusterSets.jl")

end
