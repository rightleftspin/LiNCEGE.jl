# Need to do
# 
#       prune:
#                       - Add type specifications
#
#       prune_par:
#                       - Parallelize the function properly
#                       - Add type specifications

"""
This is step two of the pipeline. In this step, the algorithm takes in
an array of clusters and a function that is used to categorize the
clusters. Using this information, the algorithm will collect clusters
that are alike by the function and combine them to output a hashmap
containing the hash of the cluster, the (potentially) modified 
representative cluster and the multiplicity of the cluster. This multiplicity
corresponds to the number of clusters that are identical to the
representative cluster under the given function.
"""

"""
Main function in step two of the pipeline. Collects clusters of the
same type and counts their corresponding multiplicities. 

Inputs: 
      clusters: array of clusters to be pruned

      pruning_fxn: function that will be applied to each cluster,
      needs to return a hash of the cluster and the cluster itself
      (the function is allowed to modify the cluster)

Output:
      Hashmap relating the hash of the cluster to the cluster itself and
      the corresponding multiplicity of the cluster.
"""
function prune(
        clusters,
        pruning_fxn
    )

    # Initialize the empty output dictionary
    cluster_mult = Dict()

    # Add function for multiplicity
    add_mult_one = (cluster, mult) -> (cluster, mult + 1)

    for cluster in clusters
        # Find the hash and rearranged cluster for each cluster
        hash, new_cluster = pruning_fxn(cluster)
        # Add this information to the output dictionary
        cluster_mult[hash] = add_mult_one(get(cluster_mult, hash, 0))
    end

    cluster_mult
end

"""
Parallel version of the main function in step two of the pipeline. Collects clusters of the
same type and counts their corresponding multiplicities. 

Inputs: 
      clusters: array of clusters to be pruned

      pruning_fxn: function that will be applied to each cluster,
      needs to return a hash of the cluster and the cluster itself
      (the function is allowed to modify the cluster)

Output:
      Hashmap relating the hash of the cluster to the cluster itself and
      the corresponding multiplicity of the cluster.
"""
function prune_par(
        clusters,
        pruning_fxn
    )

    # Initialize the empty output dictionary
    cluster_mult = Dict()

    # Add function for multiplicity
    add_mult_one = (cluster, mult) -> (cluster, mult + 1)

    for cluster in clusters
        # Find the hash and rearranged cluster for each cluster
        hash, new_cluster = pruning_fxn(cluster)
        # Add this information to the output dictionary
        cluster_mult[hash] = add_mult_one(get(cluster_mult, hash, 0))
    end

    cluster_mult
end
