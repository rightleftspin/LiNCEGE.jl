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
    center::AbstractVector{<:Integer}
    vertex_labels::AbstractVector{<:Integer}
    adj_list::AbstractVector{<:AbstractVector{<:Integer}}
    adj_matrix::AbstractMatrix{<:Integer}
    adj_matrix_weights::AbstractArray{<:Integer,3}
    coordinates::AbstractVector{<:AbstractVector{<:Real}}
    permutations::Union{<:AbstractVector, Nothing}

    function NLCELattice(
        center::AbstractVector{<:Integer},
        vertex_labels::AbstractVector{<:Integer},
        adj_list::AbstractVector{<:AbstractVector{<:Integer}},
        adj_matrix::AbstractMatrix{<:Integer},
        adj_matrix_weights::AbstractArray{<:Integer,3},
        coordinates::AbstractVector{<:AbstractVector{<:Real}},
        permutations::Union{<:AbstractVector, Nothing},
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
            coordinates,
            permutations,
        )
    end
end

# Adjacency Matrix Constructor
function NLCELattice(
    center::AbstractVector{<:Integer},
    vertex_labels::AbstractVector{<:Integer},
    adj_matrix::AbstractMatrix{<:Integer},
    adj_matrix_weights::AbstractArray{<:Integer,3},
    coordinates::AbstractVector{<:AbstractVector{<:Real}},
    permutations::Union{<:AbstractVector, Nothing},
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
        coordinates,
        permutations,
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
    basis::AbstractVector{<:AbstractVector{<:Real}},
    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    neighborhood::AbstractVector{<:Real},
    max_order::Integer;
    symmetries::Union{Nothing, <:AbstractVector} = nothing,
    basis_colors::AbstractVector{<:Integer} = repeat([1], length(basis)),
)

    coordinates, sublattice_coords, colors, centers =
        generate_coordinates(basis, primitive_vectors, max_order, basis_colors)

    permutations = nothing

    if symmetries != nothing
        if max_order > 1
            shifts = [site ./ 2 for site in eachcol(generate_primitive_lattice(primitive_vectors, 2)[1])]
            permutations = find_permutations(coordinates, symmetries, shifts)
        else
            shifts = [site ./ 2 for site in eachcol(generate_primitive_lattice(primitive_vectors, 1)[1])]
            permutations = find_permutations(coordinates, symmetries, shifts)
        end
    end

    number_vertices = length(coordinates)
    adj_matrix = zeros(Int, number_vertices, number_vertices)
    adj_matrix_weights = zeros(Int, 2, number_vertices, number_vertices)

    directions::Vector{Vector{Real}} = []

    for (index_coord, coord) in enumerate(coordinates)
        for (index_distance, distance) in enumerate(neighborhood)
            equal_distance = n -> sqrt(sum((coord - n[2]) .^ 2)) ≈ distance
            # Find all neighbors equal to the current distance but after the last distance
            for (index_neighbor, neighbor) in
                filter(equal_distance, collect(enumerate(coordinates)))
                direction = neighbor - coord
                # Check the direction of the bond, ie, along which axis
                if findfirst(≈(direction), directions) != nothing
                    adj_matrix[index_coord, index_neighbor] = 1
                    adj_matrix_weights[1, index_coord, index_neighbor] = index_distance
                    adj_matrix_weights[2, index_coord, index_neighbor] =
                        findfirst(≈(direction), directions)
                else
                    append!(directions, [direction])
                    adj_matrix[index_coord, index_neighbor] = 1
                    adj_matrix_weights[1, index_coord, index_neighbor] = index_distance
                    adj_matrix_weights[2, index_coord, index_neighbor] =
                        findfirst(≈(direction), directions)
                end
            end
        end
    end

    #NLCELattice(
    #    centers,
    #    colors,
    #    adj_matrix,
    #    adj_matrix_weights,
    #    coordinates,
    #    permutations,
    #    false,
    #    (length(neighborhood) > 1),
    #    (length(unique(basis_colors)) > 1),
    #)

    NLCELattice(
        centers,
        colors,
        adj_matrix,
        adj_matrix_weights,
        sublattice_coords,
        permutations,
        false,
        (length(neighborhood) > 1),
        (length(unique(basis_colors)) > 1),
    )

end

begin #Required functions for the pipeline
    nv(lattice::NLCELattice) = lattice.number_vertices
    center(lattice::NLCELattice) = lattice.center
    neighbors(lattice::NLCELattice, vertices::Union{Integer,AbstractArray}) =
        lattice.adj_list[vertices]
    label(lattice::NLCELattice, vertices::Union{Integer,AbstractArray}) =
        lattice.vertex_labels[vertices]
    adjacency_matrix(lattice::NLCELattice) = lattice.adj_matrix
    adjacency_matrix_weights(lattice::NLCELattice) = lattice.adj_matrix_weights

    cluster(lattice::NLCELattice, vertices::AbstractVector{<:Integer}) = NLCECluster(
        copy(vertices),
        lattice,
        adjacency_matrix(lattice)[vertices, vertices],
        adjacency_matrix_weights(lattice)[:, vertices, vertices],
    )

    get_coordinates(lattice::NLCELattice, vertices::Union{Integer,AbstractArray}) =
        lattice.coordinates[vertices]

    all_coordinates(lattice::NLCELattice) = lattice.coordinates
    permutations(lattice::NLCELattice) = lattice.permutations
end
