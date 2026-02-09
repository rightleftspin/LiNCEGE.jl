struct IsomorphicClusters{H<:IsomorphicHash,C<:IsomorphicCluster} <: AbstractClusters{H,C}
    clusters::Dict{H,C}
end

IsomorphicClusters() = IsomorphicClusters(Dict{IsomorphicHash,IsomorphicCluster}())
IsomorphicClusters{IsomorphicHash,IsomorphicCluster}() = IsomorphicClusters(Dict{IsomorphicHash,IsomorphicCluster}())

function IsomorphicClusters(translation_clusters::TranslationClusters, lattice::AbstractLattice)
    refined = IsomorphicClusters()
    rlock = ReentrantLock()

    @info "Starting parallel refinement with $(nthreads()) threads."

    for (_, translation_cluster) in translation_clusters
        cluster = IsomorphicCluster(translation_cluster, lattice)

        lock(rlock) do
            refined[ghash(cluster)] = cluster
        end
    end

    @info "Refinement complete with $(length(refined)) unique clusters."

    refined
end

_clusters(cs::IsomorphicClusters) = cs.clusters
Base.length(cs::IsomorphicClusters) = length(_clusters(cs))
Base.iterate(cs::IsomorphicClusters) = iterate(_clusters(cs))
Base.iterate(cs::IsomorphicClusters, state) = iterate(_clusters(cs), state)
Base.getindex(cs::IsomorphicClusters, ghash::IsomorphicHash) = getindex(cs.clusters, ghash)
Base.haskey(cs::IsomorphicClusters, ghash::IsomorphicHash) = haskey(cs.clusters, ghash)
function Base.setindex!(cs::IsomorphicClusters, cluster::IsomorphicCluster, ghash::IsomorphicHash)
    if ghash in cs
        cs.clusters[ghash] = merge!(cs[ghash], cluster)
    else
        cs.clusters[ghash] = cluster
    end
    cs
end
