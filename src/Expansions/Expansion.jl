"""
    Expansion(clusters, lattice, max_order)

Arbitrary expansion of clusters in the NLCE sense, contains the weights necessary to perform the NLCE summation.
"""
struct Expansion <: AbstractExpansion
        expansion_clusters::Dict{UInt,ExpansionCluster}
        order_ids::Vector{Vector{UInt}}
end

function Expansion(clusters::AbstractClusterSet, lattice::SiteExpansionLattice)
        expansion_clusters = Dict{UInt,ExpansionCluster}()
        order_ids = [Vector{UInt}() for _ in 1:max_order(lattice)]

        for cluster in clusters
                push!(order_ids[length(cluster)], cluster.ghash)
                expansion_clusters[cluster.ghash] = ExpansionCluster(cluster, clusters, lattice)
        end

        Expansion(expansion_clusters, order_ids)
end

function Expansion(clusters::AbstractClusterSet, lattice::AbstractClusterExpansionLattice)
        expansion_clusters = Dict{UInt,ExpansionCluster}()
        order_ids = [Vector{UInt}() for _ in 1:(max_order(lattice)+1)]

        # Adds Single Sites to the cluster expansion
        lv::Int = 1
        n_single_site_clusters = n_unique_sites(clusters)
        while length(order_ids[1]) < n_single_site_clusters

                lv_hash = ghash(clusters, LatticeVertices(lv))

                if !haskey(expansion_clusters, lv_hash)
                        expansion_clusters[lv_hash] = ExpansionCluster(lv_hash, n_single_site_clusters)
                        push!(order_ids[1], lv_hash)
                end

                lv += 1
        end

        for cluster in clusters
                push!(order_ids[length(cluster)+1], cluster.ghash)
                expansion_clusters[cluster.ghash] = ExpansionCluster(cluster, clusters, lattice)
        end

        Expansion(expansion_clusters, order_ids)
end

Base.getindex(e::Expansion, cluster_hash::UInt) = e.expansion_clusters[cluster_hash]
each_order(e::Expansion, max_order::Int) = e.order_ids[1:max_order]

function weights(e::Expansion, order::Int)
        cluster_hashes = e.order_ids[order]
        merge(+, [get_nlce_contribution(e.expansion_clusters[ch]) for ch in cluster_hashes]...)
end
