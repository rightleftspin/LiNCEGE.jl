# Need to do
#
#       isomorphic_pruning:
#                       - account for edge weights

"""
These are pruning functions to be used inside step 2 of the algorithm,
these functions will always take in a cluster in some form and return
a hash of the cluster and a potential permutation of the cluster. If
the function returns nothing, this means that no permutation is necessary.
"""

"""
Takes a cluster and finds the translationally invariant
form of it.

Inputs: 
      cluster: takes a cluster struct, uses the direction
      weights matrix in the cluster

Output:
      Tuple of cluster hash and nothing
"""
function symmetric_pruning(cluster::AbstractNLCECluster)

    (hash(sum(weight -> 2^weight, direction_matrix(cluster)[sortperm(vertices(cluster)), :], dims=2)), nothing)

end

"""
Finds the canonical ordering of the cluster using Nauty, this
function rearranges the cluster according to the permutation
used by Nauty

Inputs: 
      cluster: takes in a cluster struct, uses the
      edge_weighted_matrix inside it

Output:
      Tuple of cluster hash and permutation from nauty
"""
function isomorphic_pruning(cluster::AbstractNLCECluster)

    # Create an empty DenseNautyGraph
    nauty_graph = NautyGraph(edge_weighted_matrix(cluster), label(cluster, vertices(cluster)))

    # Canonize and find the corresponding permutation 
    permutation, _ = canonize!(nauty_graph)

    # Return the nauty hash and the permutation for the cluster
    # the slice is because the permutation will potentially
    # be longer than the initial graph, since edge weights
    # add extra vertices
    # Permutation goes backwards from the canonized graph to the
    # original graph
    return (ghash(nauty_graph), permutation[1:nv(cluster)])

end
