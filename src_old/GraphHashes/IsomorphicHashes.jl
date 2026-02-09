struct IsomorphicHash{H<:Unsigned} <: AbstractGraphHash{H}
    hash::H
end

function IsomorphicHash(vs::AbstractVertices, lattice::AbstractLattice; initial_ghash=nothing)

    if is_weighted(lattice)
        (h, p) = weighted_iso_hash(
            get_isomorphic_matrix(lattice, vs),
            get_labels(lattice, vs)
        )
    else
        (h, p) = unweighted_iso_hash(
            get_isomorphic_matrix(lattice, vs),
            get_labels(lattice, vs)
        )
    end

    if !isnothing(initial_ghash)
        return (IsomorphicHash(h), IsomorphicPermutation(initial_ghash, p))
    end

    (IsomorphicHash(h), EmptyPermutation())
end

IsomorphicHash{UInt}(vs::AbstractVertices, lattice::AbstractLattice; initial_ghash=nothing) = IsomorphicHash(vs, lattice, initial_ghash=initial_ghash)

"""
Finds the canonical ordering of an edge labeled cluster using Nauty, this function
returns the permutation to rearrange the cluster used by Nauty
"""
function weighted_iso_hash(weighted_adj_mat::AbstractMatrix{Int}, labels::AbstractVector{Int})

    # This is to get edge-labels to work
    nauty_labels = vcat(labels, zeros(Int64, fld(count(>(1), weighted_adj_mat), 2)))

    unweighted_adjacency_matrix = zeros(Int64, length(nauty_labels), length(nauty_labels))

    current_aux_vert = size(weighted_adj_mat, 1) + 1
    for j = 1:size(weighted_adj_mat, 2)
        for i = j:size(weighted_adj_mat, 1)

            if (@inbounds(weighted_adj_mat[i, j]) == 1)
                # Add an edge if there is already an edge
                @inbounds unweighted_adjacency_matrix[i, j] = 1
                @inbounds unweighted_adjacency_matrix[j, i] = 1
            elseif (@inbounds(weighted_adj_mat[i, j]) != 0)
                # Add an edge to the auxilary vertex here
                @inbounds unweighted_adjacency_matrix[i, current_aux_vert] = 1
                @inbounds unweighted_adjacency_matrix[j, current_aux_vert] = 1
                @inbounds unweighted_adjacency_matrix[current_aux_vert, i] = 1
                @inbounds unweighted_adjacency_matrix[current_aux_vert, j] = 1
                # Color the aux vertex the same as the edge
                @inbounds nauty_labels[current_aux_vert] = weighted_adj_mat[i, j]
                current_aux_vert += 1
            end
        end
    end

    nauty_graph = NautyGraph(unweighted_adjacency_matrix, vertex_labels=nauty_labels)
    # Canonize and find the corresponding permutation
    permutation = canonize!(nauty_graph)

    # TODO: Add cluster symmetries correctly

    # Return the nauty hash and the permutation for the cluster
    # the slice is because the permutation will potentially
    # be longer than the initial graph, since edge weights
    # add extra vertices
    # Permutation goes from the original graph to the
    # canonized graph
    return (NautyGraphs.ghash(nauty_graph), @inbounds(permutation[1:size(weighted_adj_mat, 1)]))

end

"""
Finds the canonical ordering of an edge unlabeled cluster using Nauty, this function
returns the permutation to rearrange the cluster used by Nauty
"""
function unweighted_iso_hash(unweighted_adj_mat::AbstractMatrix{Int}, labels::AbstractVector{Int})

    nauty_graph = NautyGraph(unweighted_adj_mat, vertex_labels=labels)
    # Canonize and find the corresponding permutation
    permutation = canonize!(nauty_graph)

    # TODO: Add cluster symmetries correctly

    # Return the nauty hash and the permutation for the cluster
    # the slice is because the permutation will potentially
    # be longer than the initial graph, since edge weights
    # add extra vertices
    # Permutation goes from the original graph to the
    # canonized graph
    return (NautyGraphs.ghash(nauty_graph), @inbounds(permutation[1:size(unweighted_adj_mat, 1)]))

end
