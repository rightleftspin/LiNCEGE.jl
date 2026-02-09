struct IsomorphicCluster{V<:AbstractVertices,H<:IsomorphicHash,P<:AbstractPermutation} <: AbstractCluster{V,H}
    vertices::V
    ghash::H
    permutations::Vector{P}
end

function IsomorphicCluster(translation_cluster::TranslationCluster, lattice::AbstractLattice)
    vs = vertices(translation_cluster)
    ghash, permutation = IsomorphicHash(vs, lattice)
    IsomorphicCluster(vs, ghash, EmptyPermutation[permutation])
end

vertices(cluster::IsomorphicCluster) = cluster.vertices
ghash(cluster::IsomorphicCluster) = cluster.ghash
lattice_constant(cluster::IsomorphicCluster) = length(permutations(cluster))

permutations(cluster::IsomorphicCluster) = cluster.permutations
function merge!(cluster1::IsomorphicCluster, cluster2::IsomorphicCluster)
    append!(cluster1.permutations, permutations(cluster2))
    cluster1
end
