"""
General utility functions for a variety of parts of the NLCE process, these
will generally be math heavy functions that are used often
"""

using LinearAlgebra

"""
Generates a lattice from the given basis and primitive lattice vectors, colors according to
the appropriate coloring
"""
function generate_coordinates(
    basis::AbstractVector{<:AbstractVector{<:Real}},
    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    max_order::Integer,
    basis_colors::AbstractVector{<:Integer},
)
    primitive_lattice, cartesian_lattice = generate_primitive_lattice(primitive_vectors, max_order)
    coordinates, sublattice_points, centers = add_basis(basis, primitive_lattice, cartesian_lattice)
    colors = repeat(basis_colors, size(primitive_lattice)[2])
    
    sorted_coord_perm = sortperm(coordinates)
    
    (coordinates[sorted_coord_perm], sublattice_points[sorted_coord_perm], colors[sorted_coord_perm], findfirst.(isapprox.(centers), (coordinates[sorted_coord_perm], )))
end

"""
Generates the lattice by populating the basis around each site in the primitive lattice
"""
function add_basis(basis, lattice, cartesian_lattice)

    (collect(Iterators.flatten([[site + elem for elem in basis] for site in eachcol(lattice)])), collect(Iterators.flatten([[[site..., i] for i in 1:length(basis)] for site in eachcol(cartesian_lattice)])), basis)

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
            if adj_matrix[i, j] == 1
                append!(neighbors, j)
            end
        end
        push!(adj_list, neighbors)
    end

    adj_list

end

"""
Converts adjacency matrix for a graph to an edge list"""
function adj_matrix_to_edge_list(adj_matrix::AbstractMatrix{<:Integer})

    # Find the number of vertices in the graph
    number_vertices = size(adj_matrix)[1]

    edge_list::Vector{Vector{Int64}} = []

    for i = 1:number_vertices
        for j = i:number_vertices
            if adj_matrix[i, j] != 0
                push!(edge_list, [i, j, adj_matrix[i, j]])
            end
        end
    end

    edge_list

end

"""
Finds a permutation representation of a given group of transformations
on the given coordinates. 
    
"""
function find_permutations(coordinates, group, shifts)
    permutations::Vector{Vector{Union{Int, Nothing}}} = []

    for elem in group
        max_perm = findfirst.(isapprox.([(elem * coord) + shifts[1] for coord in coordinates], atol = 1e-5), (coordinates,))

        for shift in shifts[2:end]
            perm = findfirst.(isapprox.([(elem * coord) + shift for coord in coordinates], atol = 1e-5), (coordinates,))
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

function write_to_file_coordinates(
    nlce_output::AbstractDict{AbstractNLCECluster,Vector{<:Real}},
    cluster_hashes::AbstractDict{AbstractNLCECluster, Integer},
    cluster_perms::AbstractDict{AbstractNLCECluster, Vector{<:Integer}},
    filename::AbstractString,
)

    nlce_file = open(filename, "w")

    for (cluster, mults) in nlce_output
        perm = cluster_perms[cluster]
        write(nlce_file, "$(nv(cluster)):")
        write(nlce_file, " $(cluster_hashes[cluster]):")
        for edge in edge_list(cluster)
            write(nlce_file, " $(join((findall(x -> x == edge[1], perm)[1], findall(x -> x == edge[2], perm)[1], edge[3]), ' '))")
        end
        write(nlce_file, ":")
        for coord in all_coordinates(cluster)[perm]
            write(nlce_file, " ($(join(coord, ',')))")
        end
        write(nlce_file, ": $(join(orbits(cluster)[perm], ' '))")
        write(nlce_file, ": $(join(mults, ' '))\n")
    end
    close(nlce_file)
end

function write_to_file(
    nlce_output::AbstractDict{AbstractNLCECluster,Vector{<:Real}},
    filename::AbstractString,
)

    nlce_file = open(filename, "w")

    for (cluster, mults) in nlce_output
        write(nlce_file, "$(nv(cluster)):")
        for edge in edge_list(cluster)
            write(nlce_file, " $(join(edge, ' '))")
        end
        write(nlce_file, " : $(join(mults, ' '))\n")
    end
    close(nlce_file)
end

function write_to_file_fortran(
    nlce_output::AbstractDict{AbstractNLCECluster,Vector{<:Real}},
    filename::AbstractString,
    max_order::Integer,
)

    nlce_files = [open(filename * "_$(i).txt", "w") for i = 1:max_order]
    sorted_clusters = sort(collect(keys(nlce_output)), by = nv)

    for cluster in sorted_clusters
        edges = edge_list(cluster)
        write(nlce_files[nv(cluster)], "$(length(edges))\n")
        for edge in edges
            write(nlce_files[nv(cluster)], "$(join(edge, '\t'))\n")
        end
        write(nlce_files[nv(cluster)], "\n")
    end

    for cluster in sorted_clusters
        write(nlce_files[nv(cluster)], "$(join(nlce_output[cluster], ' '))\n")
    end

    close.(nlce_files)
end
