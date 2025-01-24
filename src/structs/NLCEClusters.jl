# Need to do
#       NLCELattice:
#                       - Raise errors if adj_list and adj_matrix don't make sense
"""
"""

abstract type AbstractNLCECluster end

struct NLCECluster <: AbstractNLCECluster
    vertices::AbstractVector
    lattice::AbstractNLCELattice
    adj_matrix::AbstractMatrix{<:Integer}
    adj_matrix_weights::AbstractArray{<:Integer,3}
    orbits::Vector{<:Integer}

    # Basic Constructor
    function NLCECluster(
        vertices::AbstractVector,
        lattice::AbstractNLCELattice,
        adj_matrix::AbstractMatrix{<:Integer},
        adj_matrix_weights::AbstractArray{<:Integer,3},
    )
        return new(vertices, lattice, adj_matrix, adj_matrix_weights, Vector{Integer}())
    end
end

begin #Required functions for the pipeline
    nv(cluster::NLCECluster) = length(cluster.vertices)
    vertices(cluster::NLCECluster) = cluster.vertices
    underlying_lattice(cluster::NLCECluster) = cluster.lattice
    unit_cell_size(cluster::NLCECluster) = length(center(underlying_lattice(cluster)))
    label(cluster::NLCECluster, vertices::Union{Integer,AbstractArray}) =
        label(underlying_lattice(cluster), vertices)
    neighbors(cluster::NLCECluster, vertex::Integer) =
        filter(in(vertices(cluster)), neighbors(underlying_lattice(cluster), vertex))
    adjacency_matrix(cluster::NLCECluster) = cluster.adj_matrix
    edge_weighted_matrix(cluster::NLCECluster) = cluster.adj_matrix_weights[1, :, :]
    direction_matrix(cluster::NLCECluster) = cluster.adj_matrix_weights[2, :, :]
    cluster(cluster::NLCECluster, vertices::AbstractVector{<:Integer}) = NLCECluster(
        copy(vertices),
        underlying_lattice(cluster),
        adjacency_matrix(underlying_lattice(cluster))[vertices, vertices],
        adjacency_matrix_weights(underlying_lattice(cluster))[:, vertices, vertices],
    )
    get_coordinates(cluster::NLCECluster, vertices::Union{Integer,AbstractArray}) =
        get_coordinates(underlying_lattice(cluster), vertices)
    all_coordinates(cluster::NLCECluster) = get_coordinates(cluster, vertices(cluster))

    function permute!(cluster::NLCECluster, permutation::AbstractVector{<:Integer})

        @assert length(permutation) == nv(cluster) "Size of permutation is different than the size of the cluster"

        cluster.vertices = vertices(cluster)[permutation]

        cluster.vertices
    end

    edge_list(cluster::NLCECluster) = adj_matrix_to_edge_list(edge_weighted_matrix(cluster))

    function add_orbits!(cluster::NLCECluster, orbits::AbstractVector{<:Integer})
        if isempty(cluster.orbits)
            append!(cluster.orbits, orbits)
        end
        
        cluster.orbits
    end

    orbits(cluster::NLCECluster) = cluster.orbits

end
