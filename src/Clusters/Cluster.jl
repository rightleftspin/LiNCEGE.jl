struct Cluster <: AbstractCluster
        evs::ExpansionVertices
        lc::Float64
        ghash::UInt64
end

function Cluster(evs::ExpansionVertices, cs::AbstractClusterSet{C,H}) where {C<:AbstractCluster,H<:AbstractHasher}
        Cluster(evs, 1, ghash(cs, evs))
end

Base.length(c::Cluster) = length(c.evs)
Base.hash(c::Cluster, h::UInt) = hash(c.ghash, h)
