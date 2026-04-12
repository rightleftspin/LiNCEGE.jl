struct IsomorphicHasher{C<:Union{<:AbstractConnections,Nothing}} <: AbstractHasher
        hashing_matrix::Matrix{Int}
        connections::C
        labels::Vector{Int}
        is_weighted::Bool
end

function IsomorphicHasher(lattice::AbstractInfiniteLattice, connections::Union{<:AbstractConnections,Nothing})
        hashing_matrix = bond_matrix(lattice)
        labels = get_site_colors(lattice)
        is_weighted = length(unique(hashing_matrix)) > 2
        IsomorphicHasher(hashing_matrix, connections, labels, is_weighted)
end

IsomorphicHasher(lattice::AbstractInfiniteLattice) = IsomorphicHasher(lattice, nothing)
IsomorphicHasher(lattice::AbstractClusterExpansionLattice) = IsomorphicHasher(lattice, connections(lattice))

n_unique_sites(h::IsomorphicHasher) = length(unique(h.labels))

function ghash(h::IsomorphicHasher{StrongClusterConnections}, evs::ExpansionVertices)
        lvs = union(LatticeVertices(), h.connections[evs])
        ghash(h, lvs)
end

function ghash(h::IsomorphicHasher{WeakClusterConnections}, evs::ExpansionVertices)
        lvs, mask = h.connections[evs]
        hm = h.hashing_matrix[lvs, lvs]
        hm[mask] .= 0
        fh, _ = if h.is_weighted
                weighted_iso_hash(hm, h.labels[lvs])
        else
                unweighted_iso_hash(hm, h.labels[lvs])
        end

        fh
end

function ghash(h::IsomorphicHasher, lvs::LatticeVertices)
        fh, _ = if h.is_weighted
                weighted_iso_hash(h.hashing_matrix[lvs, lvs], h.labels[lvs])
        else
                unweighted_iso_hash(h.hashing_matrix[lvs, lvs], h.labels[lvs])
        end

        fh
end

"""
Finds the canonical ordering of an edge-labeled cluster using Nauty.
Returns `(nauty_hash, permutation)`.
"""
function weighted_iso_hash(weighted_adj_mat::AbstractMatrix{Int}, labels::AbstractVector{Int})
        # Encode edge weights by inserting auxiliary vertices (one per non-unit weighted edge).
        n_aux = fld(count(>(1), weighted_adj_mat), 2)
        nauty_labels = vcat(labels, zeros(Int64, n_aux))

        n = size(weighted_adj_mat, 1)
        unweighted_adjacency_matrix = zeros(Int64, n + n_aux, n + n_aux)

        current_aux_vert = n + 1
        for j in 1:size(weighted_adj_mat, 2), i in j:size(weighted_adj_mat, 1)
                w = @inbounds weighted_adj_mat[i, j]
                if w == 1
                        @inbounds unweighted_adjacency_matrix[i, j] = 1
                        @inbounds unweighted_adjacency_matrix[j, i] = 1
                elseif w != 0
                        @inbounds unweighted_adjacency_matrix[i, current_aux_vert] = 1
                        @inbounds unweighted_adjacency_matrix[j, current_aux_vert] = 1
                        @inbounds unweighted_adjacency_matrix[current_aux_vert, i] = 1
                        @inbounds unweighted_adjacency_matrix[current_aux_vert, j] = 1
                        @inbounds nauty_labels[current_aux_vert] = w
                        current_aux_vert += 1
                end
        end

        nauty_graph = NautyGraph(unweighted_adjacency_matrix, vertex_labels=nauty_labels)
        permutation = canonize!(nauty_graph)
        return (NautyGraphs.ghash(nauty_graph), @inbounds(permutation[1:n]))
end

"""
Finds the canonical ordering of an edge-unlabeled cluster using Nauty.
Returns `(nauty_hash, permutation)`.
"""
function unweighted_iso_hash(unweighted_adj_mat::AbstractMatrix{Int}, labels::AbstractVector{Int})
        nauty_graph = NautyGraph(unweighted_adj_mat, vertex_labels=labels)
        permutation = canonize!(nauty_graph)
        return (NautyGraphs.ghash(nauty_graph), @inbounds(permutation[1:size(unweighted_adj_mat, 1)]))
end
