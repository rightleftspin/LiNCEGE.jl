module Expansions

using JSON

import LINCEGE:
        Vertices.ExpansionVertices,
        Lattices.AbstractLattice,
        Lattices.SiteExpansionLattice,
        Lattices.neighbors,
        Clusters.AbstractCluster,
        Clusters.AbstractClusterSet,
        Clusters.ghash

abstract type AbstractExpansion end

Base.getindex(e::AbstractExpansion, cluster_id::Int, order::Int) = _NI("Base.getindex")
Base.length(e::AbstractExpansion) = _NI("Base.length")
add_array!(e::AbstractExpansion, order::Int, per_cluster::AbstractVector{Float64}) = _NI("add_array!")
order_ids(e::AbstractExpansion, order::Int) = _NI("num_clusters")
get_subclusters(e::AbstractExpansion, cluster_id::Int) = _NI("get_subclusters")

function summation!(e::AbstractExpansion, max_order::Int)
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

function write_to_file(e::AbstractExpansion, cs::AbstractClusterSet, lattice::AbstractLattice) 
        clusters = []
        for cluster in cs 
                push!(clusters, 
                        Dict(
                                "Bond List" => [],
                                "Coordinates" => [],
                                "Site Colors" => [],
                                "Multiplicities" => [],
                        )
                )
        end

        
        json_clusters = JSON.json(clusters; pretty_print=true)

end

include("util.jl")
include("SiteExpansions.jl")

end
