"""
General utility functions for a variety of parts of the NLCE process, these
will generally be math heavy functions that are used often
"""

using Distances
using LinearAlgebra
using Plots

"""
Generates the lattice by populating the basis around each site in the primitive lattice,
this adds it in terms of real space coordinates
"""
function add_basis_coords(basis, lattice)
    reduce(hcat, Iterators.flatten([[site + elem for elem in basis] for site in eachcol(lattice)]))
end

"""
Generates the lattice by populating the basis around each site in the primitive lattice,
this adds it in terms of sublattice coordinates
"""
function add_basis_sublattice(basis, unrotated_lattice)
     reduce(hcat, Iterators.flatten([[[site..., i] for i in 1:length(basis)] for site in eachcol(unrotated_lattice)]))
end

"""
Generates every point from the primitive lattice vectors, in a cube of
side length = 2 * maximum_order + 1 centered at the origin
"""
function generate_primitive_lattice(primitive_vectors, max_order)

    unrotated_coords = generate_cartesian_coordinates(length(primitive_vectors), max_order)
    # rotate and stretch standard cartesian cube coordinates into primitive lattice
    (stack(primitive_vectors) * unrotated_coords, unrotated_coords)

end

"""
Generate every integer point in a cube of TWICE the given side length, centered around
the origin twice the side length is necessary to make sure a point exists at the origin
"""
function generate_cartesian_coordinates(dimension, half_side_length)
    # Forces the lattice to have a strict center point
    diameter = 2 * half_side_length + 1
    # Total number of coordinates for the entire lattice
    max_coords = diameter ^ dimension
    coords = repeat(transpose(0:max_coords - 1), dimension)

    for dim = 0:(dimension - 1)
        coords[dim + 1, :] = div.(coords[dim + 1, :], diameter ^ dim) .% diameter .- half_side_length
    end

    coords
end

function full_adj_matrices(coords, neighbor_distances, colors, rev_connections, connections)
    num_coords = size(coords)[2]
    adj_mats = zeros(Int, 2, num_coords, num_coords)
    # Generates pairwise distances between all coordinates
    dist_matrix = pairwise(euclidean, coords, dims=2)
    directions = []

    for (i, coordi) in enumerate(eachcol(coords))
        for (j, coordj) in enumerate(eachcol(coords))
#        allowed_neighbors = collect(Iterators.flatten(connections[rev_connections[i]]))
#        for (j, coordj) in zip(allowed_neighbors, eachcol(coords)[allowed_neighbors])
        #    if (j in collect(Iterators.flatten(connections[rev_connections[i]])))
                for (index_distance, distance) in enumerate(neighbor_distances)
                    if isapprox(dist_matrix[i, j], distance)
                        adj_mats[1, i, j] = index_distance
                        direction = coordi - coordj
                        if findfirst(≈(direction), directions) != nothing
                            adj_mats[2, i, j] = findfirst(≈(direction), directions)
                        else
                            push!(directions, direction)
                            adj_mats[2, i, j] = findfirst(≈(direction), directions)
                        end
                    end
                end
         #   end
        end
    end

    (adj_mats, dist_matrix)
end

function adj_list(coords, neighbor_distances)
    num_coords = size(coords)[2]
    adj_list::Vector{Vector{Int}} = []
    # Generates pairwise distances between all coordinates
    dists = pairwise(euclidean, coords, dims=2)

    for i in 1:num_coords
        temp_adj_list::Vector{Int} = []
        for j in 1:num_coords
            for d in neighbor_distances
                if (dists[i, j] ≈ d)
                    append!(temp_adj_list, j)
                end
            end
        end
        push!(adj_list, temp_adj_list)
    end

    adj_list
end

function hashing_lattice_coords(real_space_coords, expansion_sublattice_coords, struct_per_basis, colors)
    sub_coords = []
    connection::Vector{Vector{Int}} = []
    rev_connection::Vector{Vector{Int}} = []
    all_colors = []

    for (ind_exp, exp_coord, r_coord) in zip(1:length(eachcol(expansion_sublattice_coords)), eachcol(expansion_sublattice_coords),
                                    eachcol(real_space_coords))
        temp_connection = []
        for (ind, sub_coord) in enumerate([per_basis + r_coord
                                           for per_basis in
                                               struct_per_basis[exp_coord[end]]])
            sub_coord_loc = findfirst(≈(sub_coord), sub_coords)
            if sub_coord_loc != nothing
                append!(temp_connection, sub_coord_loc)
                append!(rev_connection[sub_coord_loc], ind_exp)
            else
                push!(sub_coords, sub_coord)
                push!(all_colors, colors[exp_coord[end]][ind])
                append!(temp_connection, length(sub_coords))
                push!(rev_connection, [ind_exp])
            end
        end
        push!(connection, temp_connection)
    end

    (reduce(hcat, sub_coords), connection, rev_connection, all_colors)
end

"""
Converts adjacency list for a graph to an adjacency matrix
"""
function adj_list_to_adj_matrix(
    adj_list::AbstractVector{<:AbstractVector{<:Integer}},
    values::AbstractVector{<:AbstractVector{<:Integer}},
)

    # Find the number of vertices in the graph
    number_vertices = length(adj_list)
    # Initialize the empty adjacency matrix
    adj_matrix::Matrix{Int64} = zeros(number_vertices, number_vertices)

    # Loop over all the vertices in the adjacency list
    for (vertex, neighbors) in enumerate(adj_list)
        # Loop over all their neighbors
        for (neighbor_index, neighbor) in enumerate(neighbors)
            # Add the corresponding value to the adjacency matrix
            adj_matrix[vertex, neighbor] = values[vertex][neighbor_index]
        end
    end

    adj_matrix
end

function adj_list_to_adj_matrix(adj_list::AbstractVector{<:AbstractVector{<:Integer}})

    # Wrapper function for an adjacency list without weights. This function
    # will instead add a weight of 1 for every edge
    values = deepcopy(adj_list)
    for i = 1:length(values)
        values[i] .= 1
    end

    adj_list_to_adj_matrix(adj_list, values)
end

"""
Converts adjacency matrix for a graph to an adjacency list
"""
function adj_matrix_to_adj_list(adj_matrix::AbstractMatrix{<:Integer})

    # Find the number of vertices in the graph
    number_vertices = size(adj_matrix)[1]

    adj_list::Vector{Vector{Int64}} = []

    for i = 1:number_vertices
        neighbors = []
        for j = 1:number_vertices
            if i != j
                if adj_matrix[i, j] != 0
                    append!(neighbors, j)
                end
            end
        end
        push!(adj_list, neighbors)
    end

    adj_list

end

"""
Converts adjacency matrix for a graph to an edge list"""
function adj_matrix_to_edge_list(adj_matrix::AbstractMatrix{<:Real})

    # Find the number of vertices in the graph
    number_vertices = size(adj_matrix)[1]

    edge_list::Vector{Vector{Int64}} = []

    for i = 1:number_vertices
        for j = i:number_vertices
            if i != j
                if adj_matrix[i, j] != 0
                    push!(edge_list, [i, j, adj_matrix[i, j]])
                end
            end
        end
    end

    edge_list

end

# TODO: This might need to be fixed? It seems right, I cant figure it out rn
function reindex_adj_list(underlying_adj_list::AbstractVector{<:AbstractVector{<:Integer}}, verts::AbstractVector{<:Integer})
    if length(verts) == 1
        adj_list::Vector{Vector{Int}} = [[]]
        adj_list
    else
        adj_matrix_to_adj_list(adj_list_to_adj_matrix(underlying_adj_list)[verts, verts])
    end
end

# TODO: Fix this omg it is so slow, try pairwise distance
"""
Finds a permutation representation of a given group of transformations
on the given coordinates. 
    
"""
function find_permutations(coordinates, group, shifts)
    permutations::Vector{Vector{Union{Int, Nothing}}} = []

    for elem in group
        transformed_coords = [(elem * coord) for coord in coordinates]

#        perm_lengths = zeros(length(shifts))
#
#        for shift_i = 1:length(shifts)
#            perm_lengths[shift_i] = 
#
        max_perm = findfirst.(isapprox.([(coord + shifts[1]) for coord in transformed_coords], atol = 1e-5), (coordinates,))

        # This loop needs to be optimized, very parallelizeable
        for shift in shifts[2:end]
            perm = findfirst.(isapprox.([(coord + shift) for coord in transformed_coords], atol = 1e-5), (coordinates,))
            if length(unique(perm)) > length(unique(max_perm))
                max_perm = copy(perm)
            end
        end
        push!(
            permutations,
            max_perm
        )
    end

    println(length.(unique.(permutations)))
    
    permutations
end
