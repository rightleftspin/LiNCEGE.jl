# Need to do
#       propogate_par:
#               - write this function
#
#       grow:
#               - document the function itself
#
#       _grow_from_site:
#               - document the function itself
#               - optimize the function, it really needs it
#
#       grow_par:
#               - write this function

"""
This is step three of the pipeline. In this step, the algorithm takes in an array of clusters and a pruning function. Using this information, the algorithm will grow subclusters of the clusters by treating the clusters as the 'underlying cluster' in the grow step. It will then take these subclusters and prune them according to the pruning function. Lastly, it will collect the clusters into a hashmap containing the input clusters and a list of all their 
corresponding subclusters.
"""

"""
Main function in this connection. Takes clusters and finds all possible 
subclusters of each cluster, pruning them along the way.

Inputs: 
      clusters: hashmap of cluster hashes and their related
      clusters and multiplicities

      pruning_fxn: function that will be applied to the array of 
      subclusters, needs to return a hash of the subcluster and 
      a potential permutation

Output:
      Hashmap relating the cluster to a hashmap of its subclusters
      and their corresponding multiplicities
"""
function propogate(
    pruning_fxn::Function,
    clusters::AbstractDict{Integer,Tuple{<:AbstractNLCECluster,<:Integer,<:AbstractDict}},
)

    for (hash, (cluster, mult)) in clusters
        # Grow all possible subclusters smaller than the cluster then 
        # prune all subclusters and add them to the output dictionary
        clusters[hash] = (cluster, mult, prune_cluster(pruning_fxn, grow(cluster)))
    end

    clusters
end

"""
TODO: Parallel version of this algorithm
"""
function propogate_par() end

"""
Collects clusters of the same type and counts their corresponding 
multiplicities. Need a separate method since we don't want to
keep all the subclusters in memory

Inputs: 
      clusters: array of clusters to be pruned

      pruning_fxn: function that will be applied to each cluster,
      needs to return a hash of the cluster and the permutation 
      required to convert the cluster to it's isomorphic form

Output:
      Hashmap relating the hash of the cluster to the corresponding 
      multiplicity of the cluster.
"""
function prune_cluster(pruning::Function, clusters::AbstractVector{<:AbstractNLCECluster})

    # Initialize the empty output dictionary
    cluster_mult = Dict{Integer,Integer}()

    for cluster in clusters
        hash, _ = pruning(cluster)
        # Add this information to the output dictionary
        cluster_mult[hash] = get(cluster_mult, hash, 0) + 1
    end

    cluster_mult
end

"""
Grows clusters from the given vertices of the underlying_cluster.

Inputs: 
      underlying_cluster: AbstractNLCECluster

Output:
      Array of subclusters of the input underlying cluster
"""
function grow(underlying_cluster::AbstractNLCECluster)
    out_array::Vector{AbstractNLCECluster} = Vector()

    for max_order = 1:(nv(underlying_cluster)-1)
        guarding_set::Set{Int} = Set([])
        for vertex in vertices(underlying_cluster)
            init_neighbors::Set{Int} = Set(
                collect(
                    filter(
                        neighbor -> !(neighbor in guarding_set),
                        neighbors(underlying_cluster, vertex),
                    ),
                ),
            )
            vertices = [vertex]
            _grow_from_site(
                underlying_cluster,
                max_order,
                vertices,
                init_neighbors,
                guarding_set,
                out_array,
            )
            push!(guarding_set, vertex)
        end
    end

    out_array
end

"""
Grows the subclusters from a specific site, up till specific order
and outputs them into the out_array.

Inputs: 
      underlying_cluster: Graph with coordinates as vertex labels,
      vertex colors, and edge weights

      subclusters_vertices: array of vertices that are the current
      subcluster

      neighbors: set of vertices that are neighbors to the current
      subcluster

      guarding_set: set of vertices not to visit

      out_array: array of subclusters of the underlying_cluster

Output:
      Technically, the output is has_int_leaf, but in practice, the 
      output is the out_array that gets added to.
"""
function _grow_from_site(
    underlying_cluster::AbstractNLCECluster,
    max_order,
    subcluster_vertices::AbstractVector{V},
    current_neighbors::AbstractSet{V},
    guarding_set::AbstractSet{V},
    out_array::AbstractVector{<:AbstractNLCECluster},
) where {V<:Integer}

    if length(subcluster_vertices) == max_order
        push!(out_array, cluster(underlying_cluster, subcluster_vertices))
        return true
    end

    has_int_leaf = false
    new_guarding_set = copy(guarding_set)

    while !isempty(current_neighbors)
        neighbor = pop!(current_neighbors)
        append!(subcluster_vertices, neighbor)

        new_neighbors = copy(current_neighbors)

        for vertex in neighbors(underlying_cluster, neighbor)
            if (
                !(vertex in subcluster_vertices) &
                !(vertex in new_guarding_set) &
                !(vertex in new_neighbors)
            )

                push!(new_neighbors, vertex)
            end
        end

        if _grow_from_site(
            underlying_cluster,
            max_order,
            subcluster_vertices,
            new_neighbors,
            new_guarding_set,
            out_array,
        )
            pop!(subcluster_vertices)
            has_int_leaf = true
        else
            pop!(subcluster_vertices)
            return (has_int_leaf)
        end
        push!(new_guarding_set, neighbor)
        if (nv(underlying_cluster) - length(new_guarding_set)) < max_order
            return (has_int_leaf)
        end
    end
    return (has_int_leaf)
end

"""
TODO: Parallel version of this algorithm
"""
#function grow_par(underlying_cluster, max_order::Int64, starting_vertices::Vector{Int64}) end
