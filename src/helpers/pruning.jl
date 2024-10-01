# Need to do
# 
#       symmetric_pruning:
#                       - Add type specifications
#                       - Document function
#                       - write efficient version, there is a way to
#                       do it with less allocations.
#
#       isomorphic_pruning:
#                       - Add type specifications
#                       - account for edge weights, vertex colors
#                       and complete the permutation (note for later
#                       the permutation goes backwards from the canonized
#                       graph to the original graph)

"""
These are pruning functions to be used inside step 2 of the algorithm,
these functions will always take in a cluster in some form and return
a hash of the cluster and the cluster itself. In some cases, the
function may edit the cluster, and return this edited version.
"""

"""
Takes a cluster and finds the translationally invariant
form of it. This does not rearrange the cluster.

Inputs: 
      cluster: takes the cluster with edge weights based on
      the direction of the bond. Cluster should be of the form
      {
      vertex: [edge_weight, edge_weight, ...],
      next_vertex: ...
      }

Output:
      Tuple of cluster hash and cluster
"""
function symmetric_pruning(cluster)
    sorted_vertices = sort(collect(keys(cluster)))
    vertex_type_cluster = []
    for vertex in sorted_vertices
        vertex_type = 0
        for edge_weight in cluster[vertex] 
            vertex_type += 2 ^ (edge_weight - 1)
        end
        push!(vertex_type_cluster, vertex_type)                
    end
         
    return(cluster, vertex_type_cluster)
end

"""
Finds the canonical ordering of the cluster using Nauty, this
function rearranges the cluster according to the permutation
used by Nauty

Inputs: 
      cluster: takes the cluster in the form of a weighted
      edge list and an array of vertex colors. For example,
      the graph below represents a three vertex graph in
      which edges (1,2) and (1,3) have weight 1 and edge
      (2,3) has weight 2. Vertices 1 and 2 have color 1 and
      vertex 3 has color 2.
      (
        [[1, 2, 1], [1, 3, 1], [2, 3, 2]],
        [1, 1, 2]
      )

      *** Note: The vertices need to be within the range
      1 -> Number of vertices in the graph, inclusive.

Output:
      Tuple of cluster hash and rearranged cluster
"""
function isomorphic_pruning(cluster)
    
    # Create an empty DenseNautyGraph
    nauty_graph = NautyGraph(length(cluster[2]))

    # Populate the graph with edges
    for (v1, v2, w) in cluster[1]
        add_edge!(nauty_graph, v1, v2)
    end

    # Canonize and find the corresponding permutation 
    permutation, _ = canonize!(nauty_graph)
    
    # Return the nauty hash and the permuted cluster
    return(ghash(nauty_graph), cluster)
    
end
