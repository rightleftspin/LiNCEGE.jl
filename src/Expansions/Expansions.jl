module Expansions

using JSON3

import LINCEGE:
        _NI,
        Vertices.ExpansionVertices,
        Lattices.AbstractLattice,
        Lattices.SiteExpansionLattice,
        Lattices.StrongClusterExpansionLattice,
        Lattices.neighbors,
        Lattices.bond_matrix,
        Lattices.get_coordinates,
        Lattices.get_site_colors,
        Clusters.AbstractCluster,
        Clusters.AbstractClusterSet,
        Clusters.get_single_site_subgraphs,
        Clusters.ghash

abstract type AbstractExpansion end

Base.getindex(e::AbstractExpansion, cluster_id::Int, order::Int) = _NI("Base.getindex")
Base.length(e::AbstractExpansion) = _NI("Base.length")
add_array!(e::AbstractExpansion, order::Int, per_cluster::AbstractVector{Float64}) = _NI("add_array!")
order_ids(e::AbstractExpansion, order::Int) = _NI("order_ids")
get_subclusters(e::AbstractExpansion, cluster_id::Int) = _NI("get_subclusters")

function _summation!(e::AbstractExpansion, max_order::Int)
        stack = Vector{Tuple{Int,Float64}}()
        per_cluster = zeros(Float64, length(e))
        for order in 1:max_order
                for cluster_id in order_ids(e, order)
                        updated_lattice_constant = e[cluster_id, order]
                        push!(stack, (cluster_id, -updated_lattice_constant))
                        while !isempty(stack)
                                current_cluster_id, contribution = pop!(stack)
                                for subcluster_id in get_subclusters(e, current_cluster_id)
                                        @inbounds per_cluster[subcluster_id] += contribution
                                        push!(stack, (subcluster_id, -contribution))
                                end
                        end
                end
                add_array!(e, order, per_cluster)
                fill!(per_cluster, zero(Float64))
        end
        e
end

include("util.jl")
include("SiteExpansions.jl")
include("StrongClusterExpansions.jl")

end
