struct ClusterSet{C<:AbstractCluster,H<:AbstractHasher} <: AbstractClusterSet{C,H}
        clusters::Set{C}
        hasher::H
end

function TranslationClusterSet(lattice::AbstractLattice)
        ClusterSet{Cluster,TranslationHasher}(
                Set{Cluster}(),
                TranslationHasher(lattice)
        )
end

function IsomorphicClusterSet(lattice::AbstractLattice)
        ClusterSet{Cluster,IsomorphicHasher}(
                Set{Cluster}(),
                IsomorphicHasher(lattice)
        )
end

Base.length(cs::ClusterSet) = length(cs.clusters)
Base.in(c::C, cs::ClusterSet{C,H}) where {C<:AbstractCluster,H} = c in cs.clusters
Base.iterate(cs::ClusterSet) = iterate(cs.clusters)
Base.iterate(cs::ClusterSet, state) = iterate(cs.clusters, state)
Base.push!(cs::ClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = push!(cs.clusters, c)
Base.pop!(cs::ClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = pop!(cs.clusters, c)
Base.sort(cs::ClusterSet) = sort(collect(cs.clusters), by=length)

ghash(cs::ClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = ghash(cs.hasher, c.evs)
ghash(cs::ClusterSet, evs::ExpansionVertices) = ghash(cs.hasher, evs)
