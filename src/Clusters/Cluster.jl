struct Cluster <: AbstractCluster
        vs::AbstractVertices
        lc::Float64
        ghash::UInt64
end

function Cluster(vs::AbstractVertices, cs::AbstractClusterSet{C,H}, lattice::AbstractInfiniteLattice) where {C<:AbstractCluster,H<:AbstractHasher}
        Cluster(vs, 1 / n_unique_sites(lattice), ghash(cs, vs))
end

Base.length(c::Cluster) = length(c.vs)
Base.hash(c::Cluster, h::UInt) = hash(c.ghash, h)
lattice_constant(c::Cluster) = c.lc
