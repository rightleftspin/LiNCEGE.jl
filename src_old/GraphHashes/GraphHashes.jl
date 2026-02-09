module GraphHashes

using NautyGraphs

import LINCEGE:
    _NI,
    Vertices.AbstractVertices,
    Lattices.AbstractLattice,
    Lattices.translation_form,
    Lattices.is_weighted,
    Lattices.get_isomorphic_matrix,
    Lattices.get_labels

abstract type AbstractGraphHash{H<:Unsigned} end

Base.hash(ghash::AbstractGraphHash, h::UInt) = hash(ghash.hash, h)
Base.:(==)(ghash1::AbstractGraphHash, ghash2::AbstractGraphHash) = ghash1.hash == ghash2.hash

include("VertexHashes.jl")
include("IsomorphicHashes.jl")
include("TranslationHashes.jl")
include("Permutations.jl")

export AbstractGraphHash,
    VertexHash,
    IsomorphicHash,
    TranslationHash,
    AbstractPermutation,
    IsomorphicPermutation,
    EmptyPermutation,
    ghash,
    mapping
end
