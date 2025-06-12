"""
Stores adj_matrix in 2^dir form, and labels for each site along the diagonal
note that the labels have to start from numbers higher than the largest
direction number
"""
struct DirectionGraph{W<:AbstractVector,M<:AbstractMatrix{<:Integer}} <:
       AbstractAdjMatrixGraph
        weight_info::W
        adj_matrix::M
end

function DirectionGraph(real_space_lattice::RealSpaceLattice, tiling::Tiling)

        dist_matrix = pairwise_distance(real_space_lattice)
        adj_matrix = @MMatrix zeros(Int, size(dist_matrix))

        directions = Vector{Float64}[]

        for (i, coordi) in enumerate(coords(real_space_lattice))
                for (j, coordj) in enumerate(coords(real_space_lattice))
                        for distance in real_space_neighbors(tiling)
                                if isapprox(@inbounds(dist_matrix[i, j]), distance)
                                        direction = coordi - coordj
                                        if findfirst(≈(direction), directions) != nothing
                                                @inbounds adj_matrix[i, j] = 2^findfirst(≈(direction), directions)
                                        else
                                                push!(directions, direction)
                                                @inbounds adj_matrix[i, j] = 2^findfirst(≈(direction), directions)
                                        end
                                end
                        end
                end
        end

        @inbounds adj_matrix[diagind(adj_matrix)] =
                2 .^ (translational_labels(real_space_lattice) .+ length(directions))

        DirectionGraph(directions, adj_matrix)
end

"""
Takes a cluster and finds the translationally invariant hash of it. This assumes that each
site in the fundamental tiling unit is colored differently. In this way, the coloring of
the translational units supersedes the labeling of each site.

Please note that this function requires all the real_space_vertices to be sorted in dimensional order,
ie. x, y, z, ...
"""
function translational_hash(
        exp_v::ExpansionVertices,
        dg::DirectionGraph,
        lattice::ExpansionLattice,
        eg::ExpansionLatticeGraph,
)
        rsv = sorted_real_space_vertices(lattice, exp_v)
        # create new graph that is subgraph of dg
        new_adj_matrix = dg.adj_matrix[rsv, rsv]
        new_adj_matrix[get_mask(eg, exp_v, rsv)] .= 0

        hash(sum(new_adj_matrix, dims=2))
end
