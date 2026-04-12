"""
    AbstractCluster

Abstract base type for a single cluster — a connected subgraph of the lattice with an
associated lattice coefficient used in the NLCE summation.

Subtypes must implement:
- `Base.length(c)` — number of sites (vertices) in the cluster
- `Base.hash(c, h)` — the cluster's graph hash (`ghash`)
- `lattice_constant` — the cluster's lattice constant
"""
abstract type AbstractCluster end

"""
    AbstractClusterSet{C<:AbstractCluster, H<:AbstractHasher}

Abstract base type for a collection of unique clusters sharing a common hasher.
Clusters that are equivalent under the hasher's symmetry are merged

Subtypes must implement:
- `Base.length(cs)` — number of stored clusters
- `Base.in(c, cs)` — membership test
- `Base.iterate(cs)` / `Base.iterate(cs, state)` — iteration over clusters
- `Base.push!(cs, c)` — add a cluster
- `Base.pop!(cs, c)` — remove and return a cluster
- `ghash(cs, c)` / `ghash(cs, vs)` — delegate to the hasher
"""
abstract type AbstractClusterSet{C<:AbstractCluster,H<:AbstractHasher} end

# Cluster Methods
Base.length(c::AbstractCluster) = _NI("Base.length")
Base.hash(c::AbstractCluster, h::UInt) = _NI("Base.hash")
lattice_constant(c::AbstractCluster) = _NI("lattice_constant")

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
ghash(cs::AbstractClusterSet, vs::AbstractVertices) = _NI("ghash")
n_unique_sites(cs::AbstractClusterSet) = _NI("n_unique_sites")

include("util.jl")
include("Cluster.jl")
include("ClusterSets.jl")
