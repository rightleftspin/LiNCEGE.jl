"""
    AbstractExpansionCluster

Abstract base type for an expansion cluster. Stores the subgraph information for a cluster, along with necessary information to write the cluster to disk.

Subtypes must implement:
- `subgraphs(cluster)` — subgraph hashes of each subgraph within the cluster
- `subtract_subcluster!(cluster, subcluster)` — subtract the weight (W_p(c)) of the given subcluster from the corresponding cluster 
"""
abstract type AbstractExpansionCluster end

subgraphs(c::AbstractExpansionCluster) = _NI("subgraphs")
subtract_subcluster!(c::AbstractExpansionCluster, sc::AbstractExpansionCluster) = _NI("subtract_subcluster!")

"""
    AbstractExpansion

Abstract base type for an NLCE expansion. Stores the subgraph relationships between
clusters and the weight matrix that `summation!` populates.

Subtypes must implement:
- `Base.getindex(e, cluster_hash)` — get the cluster corresponding to the `cluster_hash`
- `each_order(e, max_order)` — vector of clusters at each order up till `max_order`
"""
abstract type AbstractExpansion end

Base.getindex(e::AbstractExpansion, cluster_hash::UInt) = _NI("getindex")
each_order(e::AbstractExpansion, max_order::Int) = _NI("each_order")

"""
    summation!(expansion, max_order)

Performs the recursive NLCE summation up till the max_order and populates the given expansion with the resultant weights
"""
function summation!(e::AbstractExpansion, max_order::Int)
        for cluster_hashes in each_order(e, max_order)
                for cluster_hash in cluster_hashes
                        cluster = e[cluster_hash]
                        for subcluster_hash in subgraphs(cluster)
                                subtract_subcluster!(cluster, e[subcluster_hash])
                        end
                end
        end
        e
end

include("util.jl")
include("ExpansionCluster.jl")
include("Expansion.jl")
