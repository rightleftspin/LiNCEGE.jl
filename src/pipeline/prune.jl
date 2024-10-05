# Need to do
# 
#       prune_par:
#                       - Parallelize the function properly

"""
This is step two of the pipeline. In this step, the algorithm takes in
an array of clusters and a function that is used to categorize the
clusters. Using this information, the algorithm will collect clusters
that are alike by the function and combine them to output a hashmap
containing the hash of the cluster, the permutation needed to modify
the cluster and the multiplicity of the cluster. This multiplicity
corresponds to the number of clusters that are identical to the
representative cluster under the given function.
"""

"""
Main function in step two of the pipeline. Collects clusters of the
same type and counts their corresponding multiplicities. 

Inputs: 
      clusters: array of clusters to be pruned

      pruning_fxn: function that will be applied to each cluster,
      needs to return a hash of the cluster and the permutation 
      required to convert the cluster to it's isomorphic form

Output:
      Hashmap relating the hash of the cluster to the cluster itself and
      the corresponding multiplicity of the cluster.
"""
function prune(pruning::Function, clusters::AbstractVector{<:AbstractNLCECluster})

    # Initialize the empty output dictionary
    cluster_mult = Dict{Integer,Tuple{<:AbstractNLCECluster,<:Integer,<:AbstractDict}}()

    # Add function for multiplicity
    add_mult_one = (cluster, mult, _) -> (cluster, mult + 1, Dict{Integer,Integer}())

    for cluster in clusters
        # Find the hash and rearranged cluster for each cluster
        hash, permutation = pruning(cluster)

        # Add this information to the output dictionary
        cluster_mult[hash] =
            add_mult_one(get(cluster_mult, hash, (cluster, 0, Dict{Integer,Integer}()))...)
    end

    cluster_mult
end

"""
TODO: Parallel version of this algorithm
"""
function prune_par() end

"""
Filters clusters that are unique  

Inputs: 
      clusters: array of clusters to be filtered

      pruning_fxn: function that will be applied to each cluster,
      needs to return a hash of the cluster and the permutation 
      required to convert the cluster to it's isomorphic form

Output:
      Array of clusters after filtering
"""
function filtering(pruning::Function, clusters::AbstractVector{<:AbstractNLCECluster})

    unique(cluster -> pruning(cluster)[1], clusters)

end
