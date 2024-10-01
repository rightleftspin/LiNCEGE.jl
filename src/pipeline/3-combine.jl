# Need to do
#       nlce_summation:
#                       - document the function itself
#                       - add type specifications
#
#       _weight:
#                       - document the function itself
#                       - add type specifications

"""
This is step three of the pipeline. In this step, the algorithm takes in
a hashmap relating cluster hashes to their multiplicities, their 
subclusters and their corresponding submultiplicities. Additionally, the
algorithm takes in the highest order to consider. Using this information, 
the algorithm recursively generates a hashmap of the clusters and their 
multiplicities. It will set the multiplicities of all clusters that are
higher than the given order to 0.
"""

"""
Main function in step three of the pipeline. Takes clusters from the given
cluster dictionary and finds their final multiplicities

Inputs: 
      clusters: Hashmap containing cluster hashes and their corresponding
      clusters, multiplicites, and subcluster information. Subcluster
      information is a hashmap of the subcluster hash and its associated
      multiplicity with reference to the cluster.

      max_order: Integer that details the maximum order that the summation
      will go to

Output:
      Hashmap of cluster hashes and their corresponding final multiplicity,
      zero for a cluster that is above the given order.
"""
function nlce_summation(clusters, order)
    cluster_weights = Dict()

    for (cluster_hash, (cluster, cluster_mult, _)) in clusters
        if nv(cluster) > order 
            cluster_weights[cluster_hash] = 0
        else
            weights = _weight(clusters, cluster_hash)
            map!(subcluster_weight -> cluster_mult * subcluster_weight, values(weights))
            cluster_weights = mergewith(+, cluster_weights, weights)
        end
    end
    
    cluster_weights
end

"""
This is the recursive function that actually calculates the NLCE weights. It is
specified in many papers, I am using page 558 eqns 5 and 6 from the paper
"A short introduction to numerical linked-cluster expansions" by Tang, Khatami, 
and Rigol. (https://arxiv.org/abs/1207.3366) 

Inputs: 
      clusters: Hashmap containing cluster hashes and their corresponding
      clusters, multiplicites, and subcluster information. Subcluster
      information is a hashmap of the subcluster hash and its associated
      multiplicity with reference to the cluster.

      cluster_hash: hash of the relevant cluster that we are finding
      the weights for

Output:
      Hashmap of cluster hashes and their corresponding multiplicity from
      the cluster specified in cluster_hash
"""
function _weight(clusters, cluster_hash)
    weight_dictionary = Dict([cluster_hash => 1])

    if nv(clusters[cluster_hash][1]) > 1
        for (subcluster_hash, (subcluster, subcluster_mult)) in clusters[cluster_hash][3]
            sub_weights = _weight(clusters, subcluster_hash)
            map!(mult -> -1 * subcluster_mult * mult, values(sub_weights))
            weight_dictionary = mergewith(+, weight_dictionary, sub_weights)
        end
    end

    weight_dictionary
end

