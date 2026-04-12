"""
    ExpansionCluster(cluster, lattice)

Cluster from a set of clusters that contains information necessary for summation! and writing to disk.
"""
struct ExpansionCluster <: AbstractExpansionCluster
        vertices::LatticeVertices
        lattice_constant::Float64
        subgraphs::Vector{UInt}
        weights::Dict{UInt,Float64}
end

function ExpansionCluster(cluster::AbstractCluster, clusters::AbstractClusterSet, lattice::SiteExpansionLattice)
        lc = lattice_constant(cluster)
        sgs = Vector{UInt}()

        for sg in get_subgraphs(cluster, lattice)
                push!(sgs, ghash(clusters, sg))
        end

        ExpansionCluster(cluster.vs, lc, sgs, Dict{UInt,Float64}(cluster.ghash => 1))
end

function ExpansionCluster(cluster::AbstractCluster, clusters::AbstractClusterSet, lattice::StrongClusterExpansionLattice)
        lc = lattice_constant(cluster)
        lvs = connections(lattice)[cluster.vs]
        sgs::Vector{UInt} = [ghash(clusters, LatticeVertices(lv)) for lv in lvs]

        for sg in get_subgraphs(cluster, lattice)
                push!(sgs, ghash(clusters, sg))
        end

        ExpansionCluster(lvs, lc, sgs, Dict{UInt,Float64}(cluster.ghash => 1))
end

function ExpansionCluster(cluster::AbstractCluster, clusters::AbstractClusterSet, lattice::WeakClusterExpansionLattice)
        lc = lattice_constant(cluster)
        lvs = just_lvs(connections(lattice), cluster.vs)
        sgs::Vector{UInt} = [ghash(clusters, LatticeVertices(lv)) for lv in lvs]

        for sg in get_subgraphs(cluster, lattice)
                push!(sgs, ghash(clusters, sg))
        end

        ExpansionCluster(lvs, lc, sgs, Dict{UInt,Float64}(cluster.ghash => 1))
end

function ExpansionCluster(lv::Int, single_site_hash::UInt, n_single_site_clusters::Int)
        lc = 1 / n_single_site_clusters

        ExpansionCluster(LatticeVertices(lv), lc, UInt[], Dict{UInt,Float64}(single_site_hash => 1))
end

subgraphs(cluster::ExpansionCluster) = cluster.subgraphs
subtract_subcluster!(cluster::ExpansionCluster, subcluster::ExpansionCluster) = merge!(+, cluster.weights, Dict{UInt,Float64}(k => -v for (k, v) in subcluster.weights))
get_nlce_contribution(cluster::ExpansionCluster) = Dict(k => cluster.lattice_constant * v for (k, v) in cluster.weights)
