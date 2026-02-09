module Clusters

using Base.Threads

import LINCEGE:
    _NI,
    Lattices.AbstractLattice

abstract type AbstractClusterID end
abstract type AbstractCluster end
abstract type AbstractClusters{C<:AbstractCluster} end

# Cluster Methods
lattice_constant(c::AbstractCluster) = _NI("lattice_constant")
order(c::AbstractCluster) = _NI("order")
get_subclusters(c::AbstractCluster) = _NI("get_subclusters")
add_subclusters!(c::AbstractCluster, cids::Vector{AbstractClusterID}) = _NI("get_subclusters")

# Cluster Collection Methods
Base.length(cs::AbstractClusters)::Int = _NI("Base.length")
Base.getindex(cs::AbstractClusters, cid::AbstractClusterID) = _NI("Base.getindex")
Base.in(cluster::C, cs::AbstractClusters{C}) where {C<:AbstractCluster} = _NI("Base.in")
Base.iterate(cs::AbstractClusters) = _NI("Base.iterate")
Base.iterate(cs::AbstractClusters, state) = _NI("Base.iterate")

get_id(cs::AbstractClusters, c::AbstractCluster) = _NI("get_id")
lattice_constant(cs::AbstractClusters, cid::AbstractClusterID) = lattice_constant(cs[cid])
order(cs::AbstractClusters, cid::AbstractClusterID) = order(cs[cid])
get_subclusters(cs::AbstractClusters, cid::AbstractClusterID) = get_subclusters(cs[cid])

Base.show(io::IO, cs::AbstractClusters) = print(io, "Cluster Collection with $(length(cs)) clusters.")

function get_order(cs::AbstractClusters, order::Int)
    filter(cid -> order(cs, cid), cs)
end

include("Subclusters/Subclusters.jl")
end
