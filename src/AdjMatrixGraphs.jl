abstract type AbstractAdjMatrixGraph end

struct BondGraph{W<:AbstractVector,M<:AbstractMatrix{<:Integer}} <: AbstractAdjMatrixGraph
        weight_info::W
        adj_matrix::M
end

"""
Finds the canonical ordering of an edge labeled cluster using Nauty, this function
returns the permutation to rearrange the cluster used by Nauty
"""
function weighted_iso_hash(g::BondGraph)

        weighted_adj_mat = adj_matrix(g)
        # This is to get edge-labels to work
        nauty_labels = vcat(
                labels(g),
                zeros(Int64, fld(count(>(1), weighted_adj_mat), 2)),
        )

        unweighted_adjacency_matrix =
                zeros(Int64, length(nauty_labels), length(nauty_labels))

        current_aux_vert = nv(cluster) + 1
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

        nauty_graph = NautyGraph(unweighted_adjacency_matrix, nauty_labels)
        # Canonize and find the corresponding permutation
        #permutation = canonize!(nauty_graph)
        permutation = canonical_permutation(nauty_graph)

        # TODO: Add cluster symmetries correctly

        # Return the nauty hash and the permutation for the cluster
        # the slice is because the permutation will potentially
        # be longer than the initial graph, since edge weights
        # add extra vertices
        # Permutation goes from the original graph to the
        # canonized graph
        return (ghash(nauty_graph), @inbounds(permutation[1:nv(g)]))

end

"""
Finds the canonical ordering of an edge unlabeled cluster using Nauty, this function
returns the permutation to rearrange the cluster used by Nauty
"""
function unweighted_iso_hash(g::BondGraph)

        nauty_graph =
                NautyGraph(adj_matrix(g), labels(g))
        # Canonize and find the corresponding permutation
        permutation = canonize!(nauty_graph)

        # TODO: Add cluster symmetries correctly

        # Return the nauty hash and the permutation for the cluster
        # the slice is because the permutation will potentially
        # be longer than the initial graph, since edge weights
        # add extra vertices
        # Permutation goes from the original graph to the
        # canonized graph
        return (ghash(nauty_graph), @inbounds(permutation[1:nv(g)])

end


function isomorphic_hash(sg::BondGraph, real_space_vertices::AbstractVector, mask::BitMatrix)

        # create new graph that is subgraph of sg
        new_adj_matrix = sg.adj_matrix[real_space_vertices, real_space_vertices]
        new_adj_matrix[mask] .= 0
        # in this case, the number of weights doesn't matter
        g = BondGraph(Bool[], new_adj_matrix)
        if nw(g) > 1
                return weighted_iso_hash(g)
        end
        
        unweighted_iso_hash(g)
end


"""
Stores adj_matrix in 2^dir form, and labels for each site along the diagonal
note that the labels have to start from numbers higher than the largest
direction number
"""
struct DirectionGraph{W<:AbstractVector,M<:AbstractMatrix{<:Integer}} <: AbstractAdjMatrixGraph
        weight_info::W
        adj_matrix::M
end

"""
Takes a cluster and finds the translationally invariant hash of it. This assumes that each
site in the fundamental tiling unit is colored differently. In this way, the coloring of
the translational units supersedes the labeling of each site.

Please note that this function requires all the real_space_vertices to be sorted in dimensional order,
ie. x, y, z, ...
"""
function translational_hash(g::DirectionGraph, real_space_vertices::AbstractVector, mask::BitMatrix)

        new_adj_matrix = g.adj_matrix[real_space_vertices, real_space_vertices]
        new_adj_matrix[mask] .= 0

        hash(sum(new_adj_matrix, dims=2))
end

struct DistanceGraph{W<:AbstractVector,M<:AbstractMatrix{<:Integer}} <: AbstractAdjMatrixGraph
        weight_info::W
        adj_matrix::M
end

distance_hash(g::DistanceGraph, real_space_vertices::AbstractVector) = hash(sum(g.adj_matrix[real_space_vertices, real_space_vertices], dims=2))

struct ExpansionLatticeGraph{M<:AbstractMatrix{<:Integer}} <: AbstractAdjMatrixGraph
        adj_matrix::M
end

function get_mask(g::ExpansionLatticeGraph, expansion_vertices::AbstractVector, real_space_vertices::AbstractVector)
        any.(in(expansion_vertices), adj_matrix(g)[real_space_vertices, real_space_vertices])
end

nv(g::AbstractAdjMatrixGraph) = size(g.adj_matrix, 1)
nw(g::AbstractAdjMatrixGraph) = length(g.weight_info)

labels(g::AbstractAdjMatrixGraph) = diag(g.adj_matrix)
adj_matrix(g::AbstractAdjMatrixGraph) = tril(g.adj_matrix, -1) + tril(g.adj_matrix, 1)

edge_list(g::AbstractAdjMatrixGraph) = adj_matrix_to_edge_list(adj_matrix(g))


Base.show(io::IO, g::AbstractAdjMatrixGraph) = print(
        io,
        "Graph with $(nv(g)) vertices, $(length(edge_list(g))) bonds, and $(nw(g)) unique weights",
)

