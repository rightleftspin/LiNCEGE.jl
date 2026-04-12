struct ClusterSet{C<:AbstractCluster,H<:AbstractHasher} <: AbstractClusterSet{C,H}
        clusters::Set{C}
        hasher::H
end

"""
    TranslationClusterSet(lattice)

ClusterSet constructor that initializes a ClusterSet with a hasher that preserves translational invariance
"""
function TranslationClusterSet(lattice::AbstractInfiniteLattice)
        ClusterSet{Cluster,TranslationHasher}(
                Set{Cluster}(),
                TranslationHasher(lattice)
        )
end

"""
    IsomorphicClusterSet(lattice)

ClusterSet constructor that initializes a ClusterSet with a hasher that preserves invariance under graph isomorphisms
"""
function IsomorphicClusterSet(lattice::AbstractLattice)
        ClusterSet{Cluster,IsomorphicHasher}(
                Set{Cluster}(),
                IsomorphicHasher(lattice)
        )
end

"""
    SymmetricClusterSet(lattice, symmetries)

ClusterSet constructor that initializes a ClusterSet with a hasher that preserves invariance under the given point group symmetries of the lattice.
"""
function SymmetricClusterSet(lattice::AbstractLattice, symmetries::Vector{Matrix{Float64}})
        ClusterSet{Cluster,SymmetricHasher}(
                Set{Cluster}(),
                SymmetricHasher(lattice, symmetries)
        )
end

SymmetricClusterSet(lattice::AbstractLattice, lattice_type::Symbol) = SymmetricClusterSet(lattice, all_lattice_symmetries[lattice_type])

Base.length(cs::ClusterSet) = length(cs.clusters)
Base.in(c::C, cs::ClusterSet{C,H}) where {C<:AbstractCluster,H} = c in cs.clusters
Base.iterate(cs::ClusterSet) = iterate(cs.clusters)
Base.iterate(cs::ClusterSet, state) = iterate(cs.clusters, state)
Base.push!(cs::ClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = push!(cs.clusters, c)
Base.pop!(cs::ClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = pop!(cs.clusters, c)
Base.sort(cs::ClusterSet) = sort(collect(cs.clusters), by=length)

ghash(cs::ClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = ghash(cs.hasher, c.vs)
ghash(cs::ClusterSet, vs::AbstractVertices) = ghash(cs.hasher, vs)
n_unique_sites(cs::ClusterSet) = n_unique_sites(cs.hasher)

Base.show(io::IO, cs::ClusterSet{C,H}) where {C,H} =
        print(io, "ClusterSet{$(nameof(H))} with $(length(cs)) clusters")
