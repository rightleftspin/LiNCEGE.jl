struct StrongClusterExpansion <: AbstractExpansion
        index_dictionary::Dict{UInt,Int}
        subgraphs::Vector{Vector{Int}}
        weights::Matrix{Float64}
        order_ids::Dict{Int,Vector{Int}}
end

function StrongClusterExpansion(clusters::AbstractClusterSet, lattice::StrongClusterExpansionLattice, max_order::Int)
        index_dictionary = Dict{UInt,Int}()
        subgraphs = Vector{Vector{Int}}()
        weights = zeros(Float64, length(clusters), max_order + 1)
        order_ids = Dict{Int,Vector{Int}}()

        single_site_clusters = get_single_site_clusters(lattice)
        for (i, cluster) in enumerate(single_site_clusters)
                if haskey(order_ids, length(cluster))
                        push!(order_ids[length(cluster)], i)
                else
                        order_ids[length(cluster)] = [i]
                end
                index_dictionary[cluster.ghash] = i
                weights[i, length(cluster)] = cluster.lc
                push!(subgraphs, Int[])
        end

        for (i, cluster) in enumerate(sort(clusters))
                if haskey(order_ids, length(cluster) + 1)
                        push!(order_ids[length(cluster)+1], i + length(single_site_clusters))
                else
                        order_ids[length(cluster)+1] = [i + length(single_site_clusters)]
                end
                index_dictionary[cluster.ghash] = i + length(single_site_clusters)
                weights[i, length(cluster)+1] = cluster.lc
                temp_subgraphs = Int[]
                for subgraph_evs in union(get_subgraphs(cluster, lattice), get_single_site_subgraphs(cluster, lattice))
                        push!(temp_subgraphs, index_dictionary[ghash(clusters, subgraph_evs)])
                end
                push!(subgraphs, temp_subgraphs)
        end

        StrongClusterExpansion(
                index_dictionary,
                subgraphs,
                weights,
                order_ids
        )
end

Base.getindex(e::StrongClusterExpansion, cluster_id::Int, order::Int) = getindex(e.weights, cluster_id, order)
Base.length(e::StrongClusterExpansion) = length(e.subgraphs)
order_ids(e::StrongClusterExpansion, order::Int) = e.order_ids[order]
get_subclusters(e::StrongClusterExpansion, cluster_id::Int) = e.subgraphs[cluster_id]
function add_array!(e::StrongClusterExpansion, order::Int, per_cluster::AbstractVector{Float64})
        @views e.weights[:, order] .+= per_cluster
end

summation!(e::StrongClusterExpansion, max_order::Int) = _summation!(e, max_order + 1)

# TODO: Fix this function
function write_to_json(e::StrongClusterExpansion, lattice::StrongClusterExpansionLattice, cs::AbstractClusterSet, filepath::String)
        all_coords = get_coordinates(lattice)
        all_colors = get_site_colors(lattice)
        adj = bond_matrix(lattice)

        clusters_data = []
        for cluster in cs
                n = length(cluster.evs)
                cluster_id = e.index_dictionary[cluster.ghash]

                coords = collect(eachcol(all_coords[:, collect(cluster.evs)]))

                colors = all_colors[cluster.evs]

                bonds = adj_mat_to_edge_list(adj[cluster.evs, cluster.evs])

                push!(clusters_data, Dict(
                        "cluster_id" => cluster_id,
                        "n_sites" => n,
                        "coordinates" => coords,
                        "site_colors" => colors,
                        "bonds" => bonds,
                        "weights" => vec(e.weights[cluster_id, :])
                ))
        end

        open(filepath, "w") do io
                JSON3.pretty(io, clusters_data)
        end
end
