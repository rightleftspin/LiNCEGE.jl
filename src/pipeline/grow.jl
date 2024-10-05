# Need to do
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
This is step one of the pipeline. In this step, the algorithm takes in
a lattice and the order that clusters should be generated till. 
Using this information, the algorithm recursively generates an array of 
clusters that are all subclusters of specified order of the lattice.
"""

"""
Main function in step one of the pipeline. Grows clusters on the lattice.

Inputs: 
      lattice: AbstractNLCELattice
      
      max_order: Integer that details the order that
      the subclusters will be generated at

Output:
      Array of subclusters of the input underlying cluster
"""
function grow(lattice::AbstractNLCELattice, max_order::Integer)
    out_array::Vector{AbstractNLCECluster} = Vector()
    guarding_set::Set{Int} = Set([])

    for vertex in center(lattice)
        init_neighbors::Set{Int} = Set(
            collect(
                filter(neighbor -> !(neighbor in guarding_set), neighbors(lattice, vertex)),
            ),
        )
        vertices = [vertex]
        _grow_from_site(
            lattice,
            max_order,
            vertices,
            init_neighbors,
            guarding_set,
            out_array,
        )
        push!(guarding_set, vertex)
    end

    out_array
end

"""
Grows the subclusters from a specific site, up till specific order
and outputs them into the out_array.

Inputs: 
      underlying_cluster: Graph with coordinates as vertex labels,
      vertex colors, and edge weights

      max_order: Integer that details the order that
      the subclusters will be generated at

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
    lattice::AbstractNLCELattice,
    max_order::Integer,
    subcluster_vertices::AbstractVector{V},
    current_neighbors::Set{V},
    guarding_set::Set{V},
    out_array::AbstractVector{<:AbstractNLCECluster},
) where {V<:Integer}

    push!(out_array, cluster(lattice, subcluster_vertices))

    if length(subcluster_vertices) == max_order
        return true
    end

    has_int_leaf = false
    new_guarding_set = copy(guarding_set)

    while !isempty(current_neighbors)
        neighbor = pop!(current_neighbors)
        append!(subcluster_vertices, neighbor)

        new_neighbors = copy(current_neighbors)

        for vertex in neighbors(lattice, neighbor)
            if (
                !(vertex in subcluster_vertices) &
                !(vertex in new_guarding_set) &
                !(vertex in new_neighbors)
            )

                push!(new_neighbors, vertex)
            end
        end

        if _grow_from_site(
            lattice,
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
        if (nv(lattice) - length(new_guarding_set)) < max_order
            return (has_int_leaf)
        end
    end
    return (has_int_leaf)
end

"""
TODO: Parallel version of this algorithm
"""
#function grow_par(underlying_cluster, max_order::Int64, starting_vertices::Vector{Int64}) end
