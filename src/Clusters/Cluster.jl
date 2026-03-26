struct Cluster <: AbstractCluster
        evs::ExpansionVertices
        lc::Float64
        ghash::UInt64
end

function Cluster(evs::ExpansionVertices, cs::AbstractClusterSet{C,H}, lattice::AbstractInfiniteLattice) where {C<:AbstractCluster,H<:AbstractHasher}
        Cluster(evs, 1 / n_unique_sites(lattice), ghash(cs, evs))
end

Base.length(c::Cluster) = length(c.evs)
Base.hash(c::Cluster, h::UInt) = hash(c.ghash, h)
get_single_site_subgraphs(c::Cluster, l::AbstractLattice) = get_single_site_subgraphs(c.evs, l)
