struct DirectionCluster{V<:ExpansionVertices,H<:Unsigned} <: AbstractCluster
        expansion_vertices::V
        neighbors::V
        ghash::H
end

function DirectionCluster(
        vertex::T,
        lattice::ExpansionLattice,
        dg::DirectionGraph,
        eg::ExpansionLatticeGraph,
) where {T<:Unsigned}
        exp_v = ExpansionVertices(vertex)

        DirectionCluster(
                exp_v,
                neighbors(lattice, vertex),
                translational_hash(
                        exp_v,
                        dg,
                        lattice,
                        eg
                ),
        )
end

function DirectionCluster(
        exp_v::ExpansionVertices,
        neighbors::ExpansionVertices,
        lattice::ExpansionLattice,
        dg::DirectionGraph,
        eg::ExpansionLatticeGraph,
)
        DirectionCluster(
                exp_v,
                neighbors,
                translational_hash(
                        exp_v,
                        dg,
                        lattice,
                        eg
                ),
        )
end

function neighbor_clusters(
        cluster::DirectionCluster,
        lattice::ExpansionLattice,
        dg::DirectionGraph,
        eg::ExpansionLatticeGraph,
)
        ns = DirectionCluster[]
        for n in cluster.neighbors
                push!(ns, neighbor_cluster(cluster, n, lattice, dg, eg))
        end
        ns
end

function neighbor_cluster(
        cluster::DirectionCluster,
        n::Vertex,
        lattice::ExpansionLattice,
        dg::DirectionGraph,
        eg::ExpansionLatticeGraph,
) where {Vertex<:Unsigned}
        DirectionCluster(
                union(cluster.expansion_vertices, n),
                union(cluster.neighbors, neighbors(lattice, n)),
                lattice,
                dg,
                eg,
        )
end
