# Need to do
#       NLCELattice:
#                       - Raise errors if adj_list and adj_matrix don't make sense
#                       for their weights
#                       - Raise errors if adj_list_weights contains weights that are zeros
#                       - Raise errors if vertex_labels contain labels that are zero

"""
This is the lattice struct that powers the entire NLCE algorithm. This struct is an immutable 
"""

abstract type AbstractNLCELattice end

struct NLCELattice{D,W,L} <: AbstractNLCELattice
    number_vertices::Integer
    center::Integer
    vertex_labels::AbstractVector{<:Integer}
    adj_list::AbstractVector{<:AbstractVector{<:Integer}}
    adj_matrix::AbstractMatrix{<:Integer}
    adj_matrix_weights::AbstractArray{<:Integer, 3}

    function NLCELattice(
        center::Integer,
        vertex_labels::AbstractVector{<:Integer},
        adj_list::AbstractVector{<:AbstractVector{<:Integer}},
        adj_matrix::AbstractMatrix{<:Integer},
        adj_matrix_weights::AbstractArray{<:Integer, 3},
        directed::Bool,
        edge_weighted::Bool,
        vertex_labeled::Bool,
    )
        # Ensuring the lattice makes sense
        n, _n = size(adj_matrix)
        m = length(adj_list)
        lab = length(vertex_labels)

        # Raising errors if it does not
        @assert n == _n "Adjacency Matrix needs to be square"
        @assert n == m "Adjacency Matrix and Adjacency list have different number of vertices"
        @assert lab == m "Wrong number of vertex labels"

        return new{directed,edge_weighted,vertex_labeled}(
            n,
            center,
            vertex_labels,
            adj_list,
            adj_matrix,
            adj_matrix_weights,
        )
    end
end

# Adjacency Matrix Constructor
function NLCELattice(
    center::Integer,
    vertex_labels::AbstractVector{<:Integer},
    adj_matrix::AbstractMatrix{<:Integer},
    adj_matrix_weights::AbstractArray{<:Integer, 3},
    directed::Bool,
    edge_weighted::Bool,
    vertex_labeled::Bool,
)

    NLCELattice(
        center,
        vertex_labels,
        adj_matrix_to_adj_list(adj_matrix),
        adj_matrix,
        adj_matrix_weights,
        directed,
        edge_weighted,
        vertex_labeled,
    )
end

"""
Constructor for NLCE lattice, takes in a basis, primitive primitive vectors and 
the maximum order and generates a lattice with the padding necessary. 
"""
function NLCELattice(
    basis::AbstractVector{<:AbstractVector{T}},
    primitive_vectors::AbstractVector{<:AbstractVector{T}},
    neighborhood::AbstractVector{T},
    max_order::Integer,
) where {T<:Real}

    max_order_padded = (2 * max_order) + 1
    center = (div(max_order_padded, 2) + 1) + ((div(max_order_padded, 2) + 1) * max_order_padded)
    coordinates = generate_coordinates(basis, primitive_vectors, max_order_padded)

    number_vertices = length(coordinates)
    adj_matrix = zeros(Int, number_vertices, number_vertices) 
    adj_matrix_weights = zeros(Int, 2, number_vertices, number_vertices) 

    directions::Vector{Vector{T}} = []

    for (index_coord, coord) in enumerate(coordinates)
        for (index_distance, distance) in enumerate(neighborhood)
            within_distance =
                n -> sqrt(sum((coord - n[2]) .^ 2)) ≈ distance
            # Find all neighbors within the current distance but after the last distance
            for (index_neighbor, neighbor) in
                filter(within_distance, collect(enumerate(coordinates)))
                direction = neighbor - coord
                # Check the direction of the bond, ie, along which axis
                if direction in directions
                    adj_matrix[index_coord, index_neighbor] = 1
                    adj_matrix_weights[1, index_coord, index_neighbor] = index_distance
                    adj_matrix_weights[2, index_coord, index_neighbor] = findfirst(==(direction), directions)
                else
                    append!(directions, [direction])
                    adj_matrix[index_coord, index_neighbor] = 1
                    adj_matrix_weights[1, index_coord, index_neighbor] = index_distance
                    adj_matrix_weights[2, index_coord, index_neighbor] = findfirst(==(direction), directions)
                end
            end
        end
    end

    NLCELattice(
                center,
                ones(Int64, number_vertices), 
                adj_matrix, 
                adj_matrix_weights, 
                false, 
                false, 
                false
               )
    
end

begin #Required functions for the pipeline
    nv(lattice::NLCELattice) = lattice.number_vertices
    center(lattice::NLCELattice) = lattice.center
    neighbors(lattice::NLCELattice, vertices::Union{Integer, AbstractArray}) = lattice.adj_list[vertices]
    label(lattice::NLCELattice, vertices::Union{Integer, AbstractArray}) = lattice.vertex_labels[vertices]
    adjacency_matrix(lattice::NLCELattice) = lattice.adj_matrix
    adjacency_matrix_weights(lattice::NLCELattice) = lattice.adj_matrix_weights

    cluster(lattice::NLCELattice, vertices::AbstractVector{<:Integer}) = NLCECluster(
            copy(vertices),
            lattice,
            adjacency_matrix(lattice)[vertices, vertices],
            adjacency_matrix_weights(lattice)[:, vertices, vertices]
        )
end
