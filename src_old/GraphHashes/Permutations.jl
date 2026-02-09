abstract type AbstractPermutation end

struct IsomorphicPermutation{H<:AbstractGraphHash} <: AbstractPermutation
    ghash::H
    mapping::Vector{Int}
end

ghash(perm::IsomorphicPermutation) = perm.ghash
mapping(perm::IsomorphicPermutation) = perm.mapping

struct EmptyPermutation <: AbstractPermutation end

export AbstractPermutation, IsomorphicPermutation, EmptyPermutation, ghash, mapping
