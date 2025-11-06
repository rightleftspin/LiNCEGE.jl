module ClusterExpansions

import DataStructures:
    Deque

import LINCEGE:
    Vertices.AbstractVertices,
    Lattices.AbstractLattice,
    Lattices.AbstractSiteExpansionLattice,
    Lattices.AbstractClusterExpansionLattice,
    Lattices.max_order,
    Lattices.n_unique_sites,
    GraphHashes.AbstractGraphHash,
    GraphHashes.IsomorphicHash,
    Clusters.AbstractCluster,
    Clusters.IsomorphicCluster,
    Clusters.Subgraph,
    Clusters.ghash,
    Clusters.lattice_constant,
    Clusters.vertices,
    Clusters.hashtype,
    Clusters.parent_hash,
    ClusterCollections.AbstractClusters,
    ClusterCollections.IsomorphicClusters,
    ClusterCollections.Subgraphs,
    ClusterCollections.get_orders,
    _NI

abstract type AbstractExpansion{H<:AbstractGraphHash,C<:AbstractCluster} end
abstract type AbstractSiteExpansion{H<:AbstractGraphHash,C<:AbstractCluster} <: AbstractExpansion{H,C} end
abstract type AbstractClusterExpansion{H<:AbstractGraphHash,C<:AbstractCluster} <: AbstractExpansion{H,C} end

Base.iterate(ce::AbstractExpansion) = _NI("iterate")
Base.iterate(ce::AbstractExpansion, state) = _NI("iterate")
Base.getindex(ce::AbstractExpansion, ghash::AbstractGraphHash, order::Int) = _NI("getindex")
Base.setindex!(ce::AbstractExpansion, lattice_constant::Real, ghash::AbstractGraphHash, order::Int) = ce.multiplicities[ghash][order] = lattice_constant
subgraphs(ce::AbstractExpansion, cluster::AbstractCluster) = _NI("subgraphs")
subgraphs(ce::AbstractExpansion, ghash::AbstractGraphHash) = _NI("subgraphs")

function Base.show(io::IO, ce::AbstractExpansion)
    println(io, "Clusters and their multiplicities:")
    for (ghash, cluster, mults) in ce
        println(io, "Cluster: ", cluster, " - Multiplicities: ", mults)
    end
end

# Try enumerating such that all the clusters are stored in a vector and not in a dictionary format
function summation!(ce::AbstractExpansion, clusters::AbstractClusters, lattice::AbstractSiteExpansionLattice)
    @info "Starting summation for orders 1 through $(max_order(lattice))."
    for (order, clusters_order) in enumerate(get_orders(clusters, lattice))
        for (ghash, cluster) in clusters_order
            updated_lattice_constant = lattice_constant(cluster) / n_unique_sites(lattice)
            ce[ghash, order] = updated_lattice_constant

            subgraphs_stack = Deque{Tuple{hashtype(cluster),Int}}()
            pushfirst!(subgraphs_stack, (ghash, 0))
            while !isempty(subgraphs_stack)
                current_ghash, current_depth = pop!(subgraphs_stack)
                for (_, subgraph) in subgraphs(ce, current_ghash)
                    phash = parent_hash(subgraph)
                    ce[phash, order] += ((-1)^(current_depth + 1)) * updated_lattice_constant
                    pushfirst!(subgraphs_stack, (phash, current_depth + 1))
                end
            end

        end
    end

    @info "summation complete."
    ce
end

include("SiteExpansions.jl")

export AbstractExpansion,
    AbstractSiteExpansion,
    AbstractClusterExpansion,
    SiteExpansion,
    summation!

end
