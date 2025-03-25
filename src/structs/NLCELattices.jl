
abstract type AbstractCluster end

"""
This is the cluster struct that powers the entire NLCE algorithm. It is designed to be very general so
it deals with many use cases.
"""
# Has underlying lattice, cluster expansion, vertex labeled, edge labeled
struct Cluster{U, C, V, E} <: AbstractCluster

    "Coordinates in the cluster, all index integers are in reference to this vector of coordinates"
    coordinates::AbstractVector{<:AbstractVector{<:Real}}
    "Centers of the cluster, where elements of the vector point to specific coordinates in the vector of coordinates above"
    center::AbstractVector{<:Integer}
    "Bonds between sites where index in the outer vector is site 1 and index in the inner vector is site 2"
    adj_list::AbstractVector{<:AbstractVector{<:Integer}}
    "Bonds between sites where the first index is a choice of representation of the cluster, second index is site 1 and third index is site 2.
    The first two representation choices are reserved for isomorphic and translational hashing respectively.
    Since there are no self loops, the diagonal of the first adjacency matrix contains the labels of each site"
    adj_matrices::AbstractArray{<:Integer,3}

    "Underlying cluster for the given cluster, or nothing if it is the underlying lattice"
    underlying_cluster::Union{AbstractCluster, Nothing}
    "Vertices that are a part of the cluster in reference to an underlying lattice"
    vertices::Union{AbstractVector{<:Integer}, Nothing}

    # Below are the fields for the cluster expansion as opposed to the site expansion
    "Bundles of coordinates where each inner vector represents a single site in the super cluster in the case of a cluster expansion"
    coordinate_bundles::Union{AbstractVector{<:AbstractVector{<:Integer}}, Nothing}
    "Bonds between sites in the super lattice. Here, the index integers are instead in reference to bundles in the coordinate bundles"
    super_adj_list::Union{AbstractVector{<:AbstractVector{<:Integer}}, Nothing}

    function Cluster(
        coordinates::AbstractVector{<:AbstractVector{<:Real}},
        center::AbstractVector{<:Integer},
        adj_list::AbstractVector{<:AbstractVector{<:Integer}},
        adj_matrices::AbstractArray{<:Integer,3},
        underlying_cluster::Union{AbstractCluster, Nothing},
        vertices::Union{AbstractVector{<:Integer}, Nothing},
        coordinate_bundles::Union{AbstractVector{<:AbstractVector{<:Integer}}, Nothing},
        super_adj_list::Union{AbstractVector{<:AbstractVector{<:Integer}}, Nothing},
        has_underlying_cluster::Bool,
        cluster_expansion::Bool,
        vertex_labeled::Bool,
        edge_labeled::Bool,
    )
        if !has_underlying_cluster
            # Ensuring the cluster makes sense, but only at the
            # highest level so that you don't do a lot of redundant
            # checks
            n, _n = size(adj_matrix)
            m = length(adj_list)
            lab = length(vertex_labels)

            # Raising errors if it does not
            @assert n == _n "Adjacency Matrix needs to be square"
            @assert n == m "Adjacency Matrix and Adjacency list have different number of vertices"
            @assert lab == m "Wrong number of vertex labels"
        end

        return new{has_underlying_cluster, cluster_expansion, vertex_labeled, edge_labeled}(
            coordinates,
            center,
            adj_list,
            adj_matrices,
            underlying_cluster,
            vertices,
            coordinate_bundles,
            super_adj_list,
        )
    end
end

"""
Basic constructor for site expansion lattice with edge weights and vertex labels
"""
function Cluster(
    coordinates::AbstractVector{<:AbstractVector{<:Real}},
    center::AbstractVector{<:Integer},
    adj_matrices::AbstractArray{<:Integer,3},
    vertex_labeled::Bool,
    edge_labeled::Bool,
)

    Cluster(
        coordinates,
        center,
        adj_matrix_to_adj_list(adj_matrices[1 : :]),
        adj_matrices,
        nothing,
        nothing,
        nothing,
        nothing,
        false,
        false,
        vertex_labeled,
        edge_labeled,
    )
end

# Constructor for cluster expansion lattice with edge weights and vertex labels

"""
Constructor for site expansion lattice, takes in a basis, primitive primitive vectors
and the maximum order and generates a cluster that represents the lattice.
"""
function Cluster(
    basis::AbstractVector{<:AbstractVector{<:Real}},
    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    neighborhood::AbstractVector{<:Real},
    max_order::Integer;
    basis_colors::AbstractVector{<:Integer} = repeat([1], length(basis)),
)

    coordinates, sublattice_coords, colors, centers =
        generate_coordinates(basis, primitive_vectors, max_order, basis_colors)

    adj_matrices = zeros(Int, 2, length(coordinates), length(coordinates))
    directions::Vector{Vector{Real}} = []

    for (index_coord, coord) in enumerate(coordinates)
        adj_matrices[1, index_coord, index_coord] = colors[index_coord]
        for (index_distance, distance) in enumerate(neighborhood)
            equal_distance = n -> sqrt(sum((coord - n[2]) .^ 2)) ≈ distance
            # Find all neighbors equal to the current distance but after the last distance
            for (index_neighbor, neighbor) in
                filter(equal_distance, collect(enumerate(coordinates)))
                direction = neighbor - coord
                # Check the direction of the bond, ie, along which axis
                if findfirst(≈(direction), directions) != nothing
                    adj_matrices[1, index_coord, index_neighbor] = index_distance
                    adj_matrices[2, index_coord, index_neighbor] =
                        findfirst(≈(direction), directions)
                else
                    append!(directions, [direction])
                    adj_matrices[1, index_coord, index_neighbor] = index_distance
                    adj_matrices[2, index_coord, index_neighbor] =
                        findfirst(≈(direction), directions)
                end
            end
        end
    end

    Cluster(
        sublattice_coords,
        centers,
        adj_matrices,
        (length(unique(basis_colors)) > 1),
        (length(neighborhood) > 1),
    )
end

"""
Takes the underlying cluster and returns a subcluster of it
"""
function Cluster(
    underlying_cluster::AbstractCluster{_, false, V, E},
    vertices::AbstractVector{<:Integer},
    )
    where {V, E}
    Cluster(
        coordinates(underlying_cluster, vertices),
        Vector(1:length(vertices)),
        [filter(in(vertices), n) for n in neighbors(underlying_cluster, vertices)],
        adj_matrices[:, vertices, vertices],
        underlying_cluster,
        vertices,
        nothing,
        nothing,
        true,
        false,
        V,
        E,
    )
end

begin # Standard Functions
    "Number of Vertices of the cluster"
    nv(cluster::Cluster) = length(cluster.coordinates)
    coordinates(cluster::Cluster, vertices::Union{Integer, AbstractArray}) =
        cluster.coordinates[vertices]
    center(cluster::Cluster) = cluster.center
    neighbors(cluster::Cluster, vertices::Union{Integer,AbstractArray}) =
        cluster.adj_list[vertices]
    adjacency_matrices(cluster::Cluster) = cluster.adj_matrices
    weighted_adjacency_matrix(cluster::Cluster) = adjacency_matrices(cluster)[1, :, :]
    direction_adjacency_matrix(cluster::Cluster) = adjacency_matrices(cluster)[2, :, :]
    label(cluster::Cluster, vertices::Union{Integer,AbstractArray}) =
        diag(weighted_adjacency_matrix(cluster))[vertices]
    underlying_cluster(cluster::Cluster{true, _, _, _}) = cluster.underlying_cluster
    vertices(cluster::Cluster{true, _, _, _}) = cluster.vertices
end
