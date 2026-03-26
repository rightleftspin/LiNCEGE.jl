struct TranslationHasher <: AbstractHasher
        hashing_matrix::Matrix{Int}
        connections::Union{<:AbstractConnections,Nothing}
end

function TranslationHasher(lattice::SiteExpansionLattice)
        pwd_matrix = pairwise_direction(get_coordinates(lattice))
        hashing_matrix, max_dir = unique_direction_indices(pwd_matrix, bond_matrix(lattice))
        diag_matrix = diagm(get_labels(lattice) .+ max_dir)
        TranslationHasher(2 .^ (hashing_matrix + diag_matrix), nothing)
end

function TranslationHasher(lattice::StrongClusterExpansionLattice)
        pwd_matrix = pairwise_direction(get_coordinates(lattice))
        hashing_matrix, max_dir = unique_direction_indices(pwd_matrix, bond_matrix(lattice))
        diag_matrix = diagm(get_labels(lattice) .+ max_dir)
        TranslationHasher(2 .^ (hashing_matrix + diag_matrix), connections(lattice))
end

function ghash(h::TranslationHasher, evs::ExpansionVertices)
        h = if isnothing(h.connections)
                hash(sum(h.hashing_matrix[evs, evs], dims=2))
        else
                lvs = union(LatticeVertices(), h.connections[evs])
                hash(sum(h.hashing_matrix[lvs, lvs], dims=2))
        end

        h
end

