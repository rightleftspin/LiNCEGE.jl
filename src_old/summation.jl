# This is the final part of the pipeline, here the goal is to take all
# the clusters and find their final multiplicities

"""
"""
function nlce_summation(clusters)
    cluster_weights = Dict()

    for (cluster_hash, (_, cluster_mult, _)) in clusters
        weights = weight(clusters, cluster_hash)
        map!(subcluster_weight -> cluster_mult * subcluster_weight, values(weights))
        cluster_weights = mergewith(+, cluster_weights, weights)
    end

    cluster_weights
end

"""
"""
function weight(clusters, cluster_hash)
    weight_dictionary = Dict([cluster_hash => 1])

    if nv(clusters[cluster_hash][1]) > 1
        for (subcluster_hash, (subcluster, subcluster_mult)) in clusters[cluster_hash][3]
            sub_weights = weight(clusters, subcluster_hash)
            map!(mult -> -1 * subcluster_mult * mult, values(sub_weights))
            weight_dictionary = mergewith(+, weight_dictionary, sub_weights)
        end
    end

    weight_dictionary
end
