struct SymmetricHasher <: AbstractHasher
        trans_hasher::TranslationHasher
        transformations::Vector{Vector{Int64}}
        connections::Union{<:AbstractConnections,Nothing}
end

function SymmetricHasher(lattice::SiteExpansionLattice, transformations::Vector{Matrix{Float64}})
        trans_hasher = TranslationHasher(2 .^ (hashing_matrix + diag_matrix), nothing)
end
