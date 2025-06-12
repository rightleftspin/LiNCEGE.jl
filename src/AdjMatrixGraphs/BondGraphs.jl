struct BondGraph{W<:AbstractVector,M<:AbstractMatrix{<:Integer}} <: AbstractAdjMatrixGraph
        weight_info::W
        adj_matrix::M
end

function BondGraph(real_space_lattice::RealSpaceLattice, tiling::Tiling)
        dist_matrix = pairwise_distance(real_space_lattice)
        adj_matrix = @MMatrix zeros(Int, size(dist_matrix))

        @inbounds adj_matrix[diagind(adj_matrix)] = labels(real_space_lattice)
        for (i, neighbor) in enumerate(real_space_neighbors(tiling))
                eq_dist = findall(isapprox(neighbor), dist_matrix)

                @inbounds adj_matrix[eq_dist] = tiling.neighbor_fn(real_space_lattice, eq_dist, i)
        end

        BondGraph(neighbors, adj_matrix)
end

"""
Finds the canonical ordering of an edge labeled cluster using Nauty, this function
returns the permutation to rearrange the cluster used by Nauty
"""
function weighted_iso_hash(g::BondGraph)

        weighted_adj_mat = adj_matrix(g)
        # This is to get edge-labels to work
        nauty_labels = vcat(labels(g), zeros(Int64, fld(count(>(1), weighted_adj_mat), 2)))

        unweighted_adjacency_matrix = @MMatrix zeros(Int64, length(nauty_labels), length(nauty_labels))

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

        nauty_graph = NautyGraph(adj_matrix(g), labels(g))
        # Canonize and find the corresponding permutation
        permutation = canonize!(nauty_graph)

        # TODO: Add cluster symmetries correctly

        # Return the nauty hash and the permutation for the cluster
        # the slice is because the permutation will potentially
        # be longer than the initial graph, since edge weights
        # add extra vertices
        # Permutation goes from the original graph to the
        # canonized graph
        return (ghash(nauty_graph), @inbounds(permutation[1:nv(g)]))

end


function isomorphic_hash(
        exp_v::ExpansionVertices,
        bg::BondGraph,
        lattice::ExpansionLattice,
        eg::ExpansionLatticeGraph,
)

        rsv = real_space_vertices(lattice, exp_v)
        # create new graph that is subgraph of sg
        new_adj_matrix = bg.adj_matrix[rsv, rsv]
        new_adj_matrix[get_mask(eg, exp_v, rsv)] .= 0
        # No need to run expensive unique function here
        g = BondGraph(Bool[], new_adj_matrix)

        if nw(bg) > 1
                return weighted_iso_hash(g)
        end

        unweighted_iso_hash(g)
end
