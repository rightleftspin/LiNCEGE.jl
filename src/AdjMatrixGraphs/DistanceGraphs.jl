struct DistanceGraph{W<:AbstractVector,M<:AbstractMatrix{<:Integer}} <:
       AbstractAdjMatrixGraph
        weight_info::W
        adj_matrix::M
end

# TODO: Add some ability to have colored bonds such that
# the symmetries of the lattice would be respected, for        
# example, the SSL has those two bonds
function DistanceGraph(real_space_lattice::RealSpaceLattice)

        dist_matrix = pairwise_distance(real_space_lattice)
        adj_matrix = @MMatrix zeros(Int, size(dist_matrix))
        neighbors = unique(dist_matrix)

        for (i, neighbor) in enumerate(neighbors)
                @inbounds adj_matrix[findall(isapprox(neighbor), dist_matrix)] .= 2^i
        end

        DistanceGraph(neighbors, adj_matrix)
end

function distance_hash(exp_v::ExpansionVertices, dg::DistanceGraph, lattice::ExpansionLattice)
        rsv = sorted_real_space_vertices(lattice, exp_v)

        hash(sum(dg.adj_matrix[rsv, rsv], dims=2))
end
