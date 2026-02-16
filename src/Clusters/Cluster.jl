struct Cluster <: AbstractCluster
    evs::ExpansionVertices
    lc::Rational
    ghash::UInt64
end

function Cluster(evs::ExpansionVertices, cs::AbstractClusterSet{C,H}) where {C<:AbstractCluster,H<:AbstractHasher}
    Cluster(evs, 1, ghash(cs, evs))
end

#function Cluster(evs::ExpansionVertices, cs::IsomorphicClusters)
#
#end

Base.length(c::Cluster) = length(c.evs)
Base.hash(c::Cluster, h::UInt) = hash(c.ghash, h)
