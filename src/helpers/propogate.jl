# Need to do
# 
#       propogate:
#                       - Add type specifications
#
#       propogate_par:
#                       - Parallelize the function properly
#                       - Add type specifications

"""
This is connection between step two and three of the pipeline. In this 
connection, the algorithm takes in an array of clusters and a pruning 
function. Using this information, the algorithm will grow subclusters
of the clusters by treating the clusters as the 'underlying cluster' 
in the grow step. It will then take these subclusters and prune them 
according to the pruning function. Lastly, it will collect the clusters 
into a hashmap containing the input clusters and a list of all their 
corresponding subclusters.
"""

"""
Main function in this connection. Takes clusters and finds all possible 
subclusters of each cluster, pruning them along the way.

Inputs: 
      clusters: array of clusters to be pruned

      pruning_fxn: function that will be applied to the array of 
      subclusters, needs to return a hash of the subcluster and 
      the subcluster itself
      (the function is allowed to modify the cluster)

Output:
      Hashmap relating the cluster to a hashmap of its subclusters
      and their corresponding multiplicities
"""
function propogate(
        clusters,
        pruning_fxn
    )

    # Initialize the empty output dictionary
    subclusters = Dict()

    for cluster in clusters
        # Grow all possible subclusters smaller than the cluster 
        # from the cluster, starting vertices are all vertices 
        all_subclusters = grow(cluster, nv(cluster) - 1, 1:nv(cluster))
        
        # Prune all subclusters and add them to the output dictionary
        subclusters[cluster] = prune(all_subclusters, pruning_fxn)
    end

    subclusters
end

"""
Parallel version of the main function in this joint. Takes clusters and finds all possible 
subclusters of each cluster, pruning them along the way.

Inputs: 
      clusters: array of clusters to be pruned

      pruning_fxn: function that will be applied to the array of 
      subclusters, needs to return a hash of the subcluster and 
      the subcluster itself
      (the function is allowed to modify the cluster)

Output:
      Hashmap relating the cluster to a hashmap of its subclusters
      and their corresponding multiplicities
"""
function propogate_par(
        clusters,
        pruning_fxn
    )

    # Initialize the empty output dictionary
    subclusters = Dict()

    for cluster in clusters
        # Grow all possible subclusters smaller than the cluster 
        # from the cluster, starting vertices are all vertices 
        all_subclusters = grow(cluster, nv(cluster) - 1, 1:nv(cluster))
        
        # Prune all subclusters and add them to the output dictionary
        subclusters[cluster] = prune(all_subclusters, pruning_fxn)
    end

    subclusters
end
