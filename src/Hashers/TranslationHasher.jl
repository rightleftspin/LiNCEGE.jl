struct TranslationHasher{C<:Union{<:AbstractConnections,Nothing}} <: AbstractHasher
        hashing_matrix::Matrix{Int}
        connections::C
end

function TranslationHasher(lattice::AbstractInfiniteLattice, connections::Union{<:AbstractConnections,Nothing})
        pwd_matrix = pairwise_direction(get_coordinates(lattice))
        hashing_matrix, max_dir = unique_direction_indices(pwd_matrix, bond_matrix(lattice))
        diag_matrix = diagm(get_labels(lattice) .+ max_dir)
        TranslationHasher(2 .^ (hashing_matrix + diag_matrix), connections)
end

TranslationHasher(lattice::AbstractInfiniteLattice) = TranslationHasher(lattice, nothing)
TranslationHasher(lattice::AbstractClusterExpansionLattice) = TranslationHasher(lattice, connections(lattice))

n_unique_sites(h::TranslationHasher) = length(unique(h.hashing_matrix[diagind(h.hashing_matrix)]))

ghash(h::TranslationHasher, lvs::LatticeVertices) = hash(sum(h.hashing_matrix[lvs, lvs], dims=2))
ghash(h::TranslationHasher{StrongClusterConnections}, evs::ExpansionVertices) = ghash(h, union(LatticeVertices(), h.connections[evs]))

function ghash(h::TranslationHasher{WeakClusterConnections}, evs::ExpansionVertices)
        lvs, mask = h.connections[evs]
        hm = h.hashing_matrix[lvs, lvs]
        hm[mask] .= 0
        hash(sum(hm, dims=2))
end
