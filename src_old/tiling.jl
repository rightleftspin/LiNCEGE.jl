"""
"""
function subcluster_alg(
    cluster,
    max_order::Int,
    subcluster_vertices::Vector{Int},
    neighbors::Set{Int},
    guarding_set::Set{Int},
    out_array,
)

    if size(subcluster_vertices)[1] == max_order
        push!(out_array, induced_subgraph(cluster, subcluster_vertices)[1])
        return true
    end

    has_int_leaf = false
    new_guarding_set = copy(guarding_set)

    while !isempty(neighbors)
        neighbor = pop!(neighbors)
        append!(subcluster_vertices, neighbor)

        new_neighbors = copy(neighbors)

        for vertex in outneighbors(cluster, neighbor)
            if (
                !(vertex in subcluster_vertices) &
                !(vertex in new_guarding_set) &
                !(vertex in new_neighbors)
            )

                push!(new_neighbors, vertex)
            end
        end

        if subcluster_alg(
            cluster,
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
        if (nv(cluster) - length(new_guarding_set)[1]) < max_order
            return (has_int_leaf)
        end
    end
    return (has_int_leaf)
end

"""
"""
function enumerate_subclusters(cluster, max_order::Int, starting_vertices::Vector{Int})

    out_array = Vector()
    guarding_set::Set{Int} = Set([])

    for vertex in starting_vertices
        neighbors::Set{Int} = Set(
            collect(
                filter(
                    neighbor -> !(neighbor in guarding_set),
                    collect(outneighbors(cluster, vertex)),
                ),
            ),
        )
        vertices = [vertex]
        subcluster_alg(cluster, max_order, vertices, neighbors, guarding_set, out_array)
        push!(guarding_set, vertex)
    end

    out_array
end
