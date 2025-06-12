begin # Growing functions to get subclusters from a cluster

    """
    The algorithm takes in a cluster and the order that subclusters should be generated till.
    Using this information, the algorithm recursively generates an array of vertices that
    are all subclusters of specified order or lower of the cluster.
    """
    function grow_lower(cluster::Cluster, start::Vector{<:Integer}, max_order::Integer)

        out_array::Vector{Vector{Int}} = Vector()
        guarding_set::Set{Int} = Set([])

        for vertex in start
            init_neighbors::Set{Int} = Set(
                collect(
                    filter(
                        neighbor -> !(neighbor in guarding_set),
                        neighbors(cluster, vertex),
                    ),
                ),
            )
            vertices = [vertex]
            _grow_from_site_lower(
                cluster,
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
    This is step one of the NLCE pipeline. The algorithm takes in
    a lattice and the order that clusters should be generated till.
    Using this information, the algorithm recursively generates an array of
    clusters that are all subclusters of specified order or lower of the lattice.
    """
    function grow_exact(cluster::Cluster, start::Vector{<:Integer}, max_order::Integer)

        out_array::Vector{Vector{Int}} = Vector()
        guarding_set::Set{Int} = Set([])

        for vertex in start
            init_neighbors::Set{Int} = Set(
                collect(
                    filter(
                        neighbor -> !(neighbor in guarding_set),
                        neighbors(cluster, vertex),
                    ),
                ),
            )
            vertices = [vertex]
            _grow_from_site_exact(
                cluster,
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
    Grows the subclusters from a specific site, up till specific order and outputs them into
    the out_array. This adds all subclusters at or below the given max_order.
    """
    function _grow_from_site_lower(
        cluster::Cluster,
        max_order::Integer,
        subcluster_vertices::AbstractVector{V},
        current_neighbors::Set{V},
        guarding_set::Set{V},
        out_array::AbstractVector{<:AbstractVector{V}},
    ) where {V<:Integer}

        push!(out_array, deepcopy(subcluster_vertices))

        if length(subcluster_vertices) == max_order
            return true
        end

        has_int_leaf = false
        new_guarding_set = copy(guarding_set)

        while !isempty(current_neighbors)
            neighbor = pop!(current_neighbors)
            append!(subcluster_vertices, neighbor)

            new_neighbors = copy(current_neighbors)

            for vertex in neighbors(cluster, neighbor)
                if (
                    !(vertex in subcluster_vertices) & !(vertex in new_guarding_set) &
                    !(vertex in new_neighbors)
                )

                    push!(new_neighbors, vertex)
                end
            end

            if _grow_from_site_lower(
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

            if (nsv(cluster) - length(new_guarding_set)) < max_order
                return (has_int_leaf)
            end
        end
        return (has_int_leaf)
    end

    """
    Grows the subclusters from a specific site, up till specific order and outputs them into
    the out_array. This adds all subclusters at the given max_order.
    """
    function _grow_from_site_exact(
        cluster::Cluster,
        max_order::Integer,
        subcluster_vertices::AbstractVector{V},
        current_neighbors::Set{V},
        guarding_set::Set{V},
        out_array::AbstractVector{<:AbstractVector{V}},
    ) where {V<:Integer}


        if length(subcluster_vertices) == max_order
            push!(out_array, deepcopy(subcluster_vertices))
            return true
        end

        has_int_leaf = false
        new_guarding_set = copy(guarding_set)

        while !isempty(current_neighbors)
            neighbor = pop!(current_neighbors)
            append!(subcluster_vertices, neighbor)

            new_neighbors = copy(current_neighbors)

            for vertex in neighbors(cluster, neighbor)
                if (
                    !(vertex in subcluster_vertices) & !(vertex in new_guarding_set) &
                    !(vertex in new_neighbors)
                )

                    push!(new_neighbors, vertex)
                end
            end

            if _grow_from_site_exact(
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

            if (nsv(cluster) - length(new_guarding_set)) < max_order
                return (has_int_leaf)
            end
        end
        return (has_int_leaf)
    end

end
