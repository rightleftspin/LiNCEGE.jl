
struct ExpansionLatticeGraph{M<:AbstractMatrix{<:Integer}} <: AbstractAdjMatrixGraph
        adj_matrix::M
end

function ExpansionLatticeGraph(real_space_lattice::RealSpaceLattice, tiling::Tiling)

        dist_matrix = pairwise_distance(real_space_lattice)
        adj_matrix = zeros(Int, size(dist_matrix))

        for ind in findall(d -> any(isapprox.(d, real_space_neighbors(tiling))), dist_matrix)

                @inbounds adj_matrix[ind] = matching_exp_vertex(real_space_lattice, ind[1], ind[2])
        end

        ExpansionLatticeGraph(adj_matrix)
end

function get_mask(
        exp_v::ExpansionVertices,
        eg::ExpansionLatticeGraph,
        lattice::ExpansionLattice
)
        rsv = real_space_vertices(lattice, exp_v)

        mask = any.(
                !in(exp_v),
                adj_matrix(eg)[rsv, rsv],
        )
        # Set diagonal to 0 so that it doesn't cancel out later
        mask[diagind(mask)] .= 0
        mask
end
