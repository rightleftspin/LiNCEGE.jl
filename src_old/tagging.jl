using Nauty
using Graphs, MetaGraphsNext

"""
"""
function symmetric_tagging(cluster::AbstractGraph{T}) where {T}
    sorted_vertices = sort(collect(labels(cluster)))
    vertex_type_cluster = []
    for (index, vertex) in enumerate(sorted_vertices)
        vertex_type = 0
        for neighbor in neighbor_labels(cluster, vertex)
            vertex_type += 2^(cluster[vertex, neighbor][2] - 1)
        end
        push!(vertex_type_cluster, vertex_type)
    end

    return (hash(vertex_type_cluster))
end

"""
"""
function isomorphic_tagging(cluster::AbstractGraph{T}) where {T}
    hash(Nauty.baked_canonical_form(SimpleGraphFromIterator(edges(cluster))).canong)
end

"""
"""
function counting(clusters, tagging_function)
    filtered_clusters = Dict()
    for cluster in clusters
        tag = tagging_function(cluster)
        cluster_entry = get(filtered_clusters, tag, (cluster, 0))
        filtered_clusters[tag] = (cluster_entry[1], cluster_entry[2] + 1)
    end

    # looks like {tag: (cluster, multiplicity)}
    filtered_clusters
end

"""
"""
function reducing(clusters, tagging_function)
    filtered_clusters = counting(clusters, tagging_function)
    reduced_clusters = Vector()
    for (tag, (cluster, multiplicity)) in filtered_clusters
        push!(reduced_clusters, cluster)
    end

    reduced_clusters
end
