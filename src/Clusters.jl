
abstract type AbstractCluster end
# Thought: Change the default hasher on a cluster to be the translational
# hash of the cluster. This means that you will not have to worry about
# having to prune a list, you just need to add all the clusters to a set
# of clusters during the grow step. This removes the need to filter
# cluster all togehter, and you can just group clusters instead.

"""
This is the cluster struct that powers the entire NLCE algorithm. It is designed to be very general so
it deals with many use cases.
"""
# Has underlying lattice, cluster expansion, vertex labeled, edge labeled
struct Cluster{U, V, E} <: AbstractCluster

    "Coordinates in the cluster, all index integers are in reference to this vector of coordinates"
    coordinates::AbstractMatrix{<:Real}
    "Start point of the cluster for growing subclusters, where elements of the vector point to specific coordinates in the vector of coordinates above
    for clustered expansions, this will refer instead to sites in the coordinate bundles"
    start::AbstractVector{<:Integer}
    "Bonds between sites where the first index is a choice of representation of the cluster, second index is site 1 and third index is site 2.
    The first two representation choices are reserved for isomorphic and translational hashing respectively.
    Since there are no self loops, the diagonal of the first adjacency matrix contains the labels of each site"
    adj_matrices::AbstractArray{<:Integer,3}
    dist_matrix::AbstractMatrix{<:Real}

    "Underlying cluster for the given cluster, or nothing if it is the underlying lattice"
    underlying_cluster::Union{AbstractCluster, Nothing}
    "Vertices that are a part of the cluster in reference to an underlying lattice"
    underlying_vertices::Union{AbstractVector{<:Integer}, Nothing}

    # Below are the fields for the cluster expansion as opposed to the site expansion
    "Bundles of coordinates where each inner vector represents a single site in the super cluster in the case of a cluster expansion"
    coordinate_bundles::Union{AbstractVector{<:AbstractVector{<:Integer}}, Nothing}
    "Bonds between sites in the super lattice. Here, the index integers are instead in reference to bundles in the coordinate bundles"
    super_adj_list::Union{AbstractVector{<:AbstractVector{<:Integer}}, Nothing}

    function Cluster(
        coordinates::AbstractMatrix{<:Real},
        start::AbstractVector{<:Integer},
        adj_matrices::AbstractArray{<:Integer,3},
        dist_matrix::AbstractMatrix{<:Real},
        underlying_cluster::Union{AbstractCluster, Nothing},
        underlying_vertices::Union{AbstractVector{<:Integer}, Nothing},
        coordinate_bundles::Union{AbstractVector{<:AbstractVector{<:Integer}}, Nothing},
        super_adj_list::Union{AbstractVector{<:AbstractVector{<:Integer}}, Nothing},
        has_underlying_cluster::Bool,
        vertex_labeled::Bool,
        edge_labeled::Bool,
    )
        #if !has_underlying_cluster
        #    # Ensuring the cluster makes sense, but only at the
        #    # highest level so that you don't do a lot of redundant
        #    # checks
        #    n, _n = size(adj_matrices[1, :, :])
        #    m = length(adj_list)

        #    # Raising errors if it does not
        #    @assert n == _n "Adjacency Matrix needs to be square"
        #    @assert n == m "Adjacency Matrix and Adjacency list have different number of vertices"
        #end

        return new{has_underlying_cluster, vertex_labeled, edge_labeled}(
            coordinates,
            start,
            adj_matrices,
            dist_matrix,
            underlying_cluster,
            underlying_vertices,
            coordinate_bundles,
            super_adj_list,
        )
    end
end

"""
Basic constructor for either expansion lattice with edge weights and vertex labels
"""
function Cluster(
    coordinates::AbstractMatrix{<:Real},
    start::AbstractVector{<:Integer},
    adj_matrices::AbstractArray{<:Integer,3},
    dist_matrix::AbstractMatrix{<:Real},
    coordinate_bundles::Union{AbstractVector{<:AbstractVector{<:Integer}}, Nothing},
    super_adj_list::Union{AbstractVector{<:AbstractVector{<:Integer}}, Nothing},
    vertex_labeled::Bool,
    edge_labeled::Bool,
)

    Cluster(
        coordinates,
        start,
        adj_matrices,
        dist_matrix,
        nothing,
        nothing,
        coordinate_bundles,
        super_adj_list,
        false,
        vertex_labeled,
        edge_labeled,
    )
end
#"""
#Constructor for a site expansion lattice, wrapper for the general cluster function
#"""
#function Cluster(
#    basis::AbstractVector{<:AbstractVector{<:Real}},
#    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
#    neighbors::AbstractVector{<:Real},
#    max_order::Integer;
#    basis_colors::AbstractVector{<:AbstractVector{<:Integer}} = repeat([1], length(basis)),
#)
#    # origin in general dimension
#    exp_basis = repeat([0], length(basis[1]))
#    #
#
#    Cluster(exp_basis, [basis],  )
#end

"""
Constructor for a lattice, generates a corresponding expansion lattice
from the specified basis and primitive vectors, where each site in the basis is the
corresponding cluster in struct_per_basis.
"""
function Cluster(
    expansion_basis::AbstractVector{<:AbstractVector{<:Real}},
    struct_per_basis::AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}},
    expansion_primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    expansion_neighbors::AbstractVector{<:Real},
    neighbors::AbstractVector{<:Real},
    max_order::Integer;
    basis_colors::AbstractVector{<:AbstractVector{<:Integer}} = [repeat([1], length(st)) for st in struct_per_basis],
)
    lattice, unrotated_lattice = generate_primitive_lattice(expansion_primitive_vectors, max_order)

    exp_coords = add_basis_coords(expansion_basis, lattice)
    exp_sublattice_coords = add_basis_sublattice(expansion_basis, unrotated_lattice)

    coords, connections, colors = hashing_lattice_coords(exp_coords,
                                                         exp_sublattice_coords,
                                                         struct_per_basis,
                                                         basis_colors)

    exp_adj_list = adj_list(coords, expansion_neighbors)
    adj_matrices, dist_matrix = full_adj_matrices(coords, neighbors, colors)
    start_points = findfirst.(isapprox.(expansion_basis), (eachcol(exp_coords), ))

    Cluster(
        coords,
        start_points,
        adj_matrices,
        dist_matrix,
        connections,
        exp_adj_list,
        (length(unique(Iterators.flatten(basis_colors))) > 1),
        (length(neighbors) > 1),
    )
end

"""
Takes the underlying cluster and returns a subcluster of it. This function is for
the cluster expansion
"""
function Cluster(
    underlying_cluster::Cluster{<:Any, V, E},
    underlying_super_vertices::AbstractVector{<:Integer},
    ) where {V, E}
    underlying_vertices = all_vertices(underlying_cluster, underlying_super_vertices)
    # Need to reindex the coordinate bundle here, this is more complex than
    # expected
    Cluster(
        coordinates(underlying_cluster, underlying_vertices),
        # This points to the cluster bundle, not the regular adjacency list!
        Vector(1:length(underlying_super_vertices)),
        adjacency_matrices(underlying_cluster)[:, underlying_vertices, underlying_vertices],
        distance_matrix(underlying_cluster)[underlying_vertices, underlying_vertices],
        underlying_cluster,
        underlying_vertices,
        reindexed_coordinate_bundle(underlying_cluster, underlying_super_vertices, Dict((underlying_vertices[i] => i for i=1:length(underlying_vertices)))),
        super_adj_list(underlying_cluster, underlying_super_vertices),
        true,
        V,
        E,
    )
end

begin # Standard Access functions
    nv(cluster::Cluster) = size(cluster.coordinates)[2]
    nv_sv(cluster::Cluster) = length(cluster.coordinate_bundles)
    nv_weights(cluster::Cluster) = (nv(cluster) == 1) ? 1 : (length(cluster.coordinate_bundles) + 1)
    dimensions(cluster::Cluster) = size(cluster.coordinates)[1]
    has_underlying_cluster(cluster::Cluster{U, <:Any, <:Any}) where {U} = U
    vertex_labeled(cluster::Cluster{<:Any, V, <:Any}) where {V} = V
    edge_labeled(cluster::Cluster{<:Any, <:Any, E}) where {E} = E
    vertices(cluster::Cluster) = Vector(1:nv(cluster))
    coordinates(cluster::Cluster, vertices::Union{Integer, AbstractVector}) = cluster.coordinates[:, vertices]
    all_coordinates(cluster::Cluster) = cluster.coordinates
    start(cluster::Cluster) = cluster.start

    neighbors(cluster::Cluster, vertices::Union{Integer,AbstractVector}) = cluster.super_adj_list[vertices]

    edge_list(cluster::Cluster) = adj_matrix_to_edge_list(weighted_adjacency_matrix(cluster))

    adjacency_matrices(cluster::Cluster) = cluster.adj_matrices
    weighted_adjacency_matrix(cluster::Cluster) = adjacency_matrices(cluster)[1, :, :] - diagm(all_labels(cluster))
    direction_adjacency_matrix(cluster::Cluster) = adjacency_matrices(cluster)[2, :, :]
    distance_matrix(cluster::Cluster) = cluster.dist_matrix

    label(cluster::Cluster, vertices::Union{Integer,AbstractVector}) = Int.(diag(weighted_adjacency_matrix(cluster))[vertices])
    all_labels(cluster::Cluster) = Int.(diag(adjacency_matrices(cluster)[1, :, :]))

    underlying_cluster(cluster::Cluster{true, <:Any, <:Any}) = cluster.underlying_cluster
    underlying_vertices(cluster::Cluster{true, <:Any, <:Any}) = cluster.underlying_vertices

    all_vertices(cluster::Cluster, super_verts::Union{Integer, AbstractVector}) = unique(vcat(coordinate_bundle(cluster, super_verts)...))
    coordinate_bundle(cluster::Cluster, super_verts::Union{Integer, AbstractVector}) = cluster.coordinate_bundles[super_verts]
    reindexed_coordinate_bundle(cluster::Cluster, super_verts::Union{Integer, AbstractVector}, reindexed_underlying_verts::AbstractDict) = [[reindexed_underlying_verts[vert] for vert in con] for con in cluster.coordinate_bundles[super_verts]]
    super_adj_list(cluster::Cluster, super_verts::Union{Integer, AbstractVector}) = reindex_adj_list(cluster.super_adj_list, super_verts)

    Base.show(io::IO, cluster::Cluster) = print(io, "Cluster with $(nv(cluster)) vertices and $(length(edge_list(cluster))) bonds. Super lattice contains $(length(cluster.super_adj_list)) super vertices")

    # Sets default hashing of a cluster to be the translationally in[1]variant hash
    Base.hash(cluster::Cluster, h::UInt) = hash(translational_pruning(cluster), h)
    Base.isequal(cluster1::Cluster, cluster2::Cluster) = (translational_pruning(cluster1) == translational_pruning(cluster2))
end

begin # Hashing Functions

    """
    Takes a vertex colored cluster and finds the
    translationally invariant hash of it.

    Inputs:
          cluster: takes a cluster struct, uses the direction weights matrix in the cluster

    Output:
          Tuple of cluster hash and nothing
    """
    function translational_pruning(cluster::Cluster{true, true, <:Any})
        perm = sortperm(eachcol(all_coordinates(cluster)))
        form = vec(
            sum(
                weight -> 2^weight,
                direction_adjacency_matrix(cluster)[perm, :],
                dims = 2,
            )
        )
        (hash(form + (1 .// (all_labels(cluster)[perm] .+ 1))), nothing)
    end

    """
    Takes an uncolored cluster and finds the
    translationally invariant hash of it.

    Inputs:
          cluster: takes a cluster struct, uses the direction weights matrix in the cluster

    Output:
          Tuple of cluster hash and nothing
    """
    function translational_pruning(cluster::Cluster{true, false, <:Any})
        (hash(vec(
            sum(
                weight -> 2^weight,
                direction_adjacency_matrix(cluster)[sortperm(eachcol(all_coordinates(cluster))), :],
                dims = 2,
            )
        )), nothing)
    end

    """
    Finds the canonical ordering of an edge-labeled, vertex colored cluster
    using Nauty, this function rearranges the cluster according to the permutation
    used by Nauty

    Inputs:
          cluster: takes in a cluster struct, uses the
                weighted_adjacency_matrix inside it

    Output:
          Tuple of cluster hash and permutation from nauty
    """
    function isomorphic_pruning(cluster::Cluster{true, <:Any, true})

        # This is to get edge-labels to work
        nauty_labels = vcat(all_labels(cluster),
                            zeros(Int64, fld(count(>(1), weighted_adjacency_matrix(cluster)), 2)))

        unweighted_adjacency_matrix = zeros(Int64, length(nauty_labels), length(nauty_labels))

        current_aux_vert = nv(cluster) + 1
        for j = 1:size(weighted_adjacency_matrix(cluster),2)
            for i = j:size(weighted_adjacency_matrix(cluster),1)

                if (weighted_adjacency_matrix(cluster)[i, j] == 1)
                    # Add an edge if there is already an edge
                    unweighted_adjacency_matrix[i, j] = 1
                    unweighted_adjacency_matrix[j, i] = 1
                elseif (weighted_adjacency_matrix(cluster)[i, j] !=0)
                    # Add an edge to the auxilary vertex here
                    unweighted_adjacency_matrix[i, current_aux_vert] = 1
                    unweighted_adjacency_matrix[j, current_aux_vert] = 1
                    unweighted_adjacency_matrix[current_aux_vert, i] = 1
                    unweighted_adjacency_matrix[current_aux_vert, j] = 1
                    # Color the aux vertex the same as the edge
                    nauty_labels[current_aux_vert] = weighted_adjacency_matrix(cluster)[i, j]
                    current_aux_vert += 1
                end
            end
        end

        nauty_graph =
            NautyGraph(unweighted_adjacency_matrix, nauty_labels)
        # Canonize and find the corresponding permutation
        #permutation = canonize!(nauty_graph)
        permutation = canonical_permutation(nauty_graph)

        # TODO: Add cluster symmetries correctly

        # Return the nauty hash and the permutation for the cluster
        # the slice is because the permutation will potentially
        # be longer than the initial graph, since edge weights
        # add extra vertices
        # Permutation goes from the original graph to the
        # canonized graph
        return (ghash(nauty_graph), permutation[1:nv(cluster)])

    end

    """
    Finds the canonical ordering of a vertex colored cluster
    using Nauty, this function rearranges the cluster according to the permutation
    used by Nauty

    Inputs:
          cluster: takes in a cluster struct, uses the
                weighted_adjacency_matrix inside it

    Output:
          Tuple of cluster hash and permutation from nauty
    """
    function isomorphic_pruning(cluster::Cluster{true, <:Any, false})

        nauty_graph =
            NautyGraph(weighted_adjacency_matrix(cluster), all_labels(cluster))
        # Canonize and find the corresponding permutation
        permutation = canonize!(nauty_graph)

        # TODO: Add cluster symmetries correctly

        # Return the nauty hash and the permutation for the cluster
        # the slice is because the permutation will potentially
        # be longer than the initial graph, since edge weights
        # add extra vertices
        # Permutation goes from the original graph to the
        # canonized graph
        return (ghash(nauty_graph), permutation[1:nv(cluster)])

    end

   # TODO: Need to fix the symmetric pruning function
   # function symmetric_pruning(cluster_prehash::AbstractNLCECluster)
   #     (
   #         hash(
   #             sort(([
   #                 translational_form(
   #                     cluster(underlying_lattice(cluster_prehash), Vector{Integer}(perm[vertices(cluster_prehash)])),
   #                 ) for perm in permutations(underlying_lattice(cluster_prehash))
   #                     ])),
   #         ),
   #         nothing,
   #     )
   # end

end

begin # Growing functions to get subclusters from a cluster
# TODO: Document all these functions, then optimize or parallelize them in some way

"""
This is step one of the NLCE pipeline. In this step, the algorithm takes in
a lattice and the order that clusters should be generated till.
Using this information, the algorithm recursively generates an array of
clusters that are all subclusters of specified order of the lattice.
"""
function grow(underlying_cluster::Cluster, max_order::Integer)
    out_array::Vector{AbstractCluster} = Vector()
    guarding_set::Set{Int} = Set([])

    for vertex in start(underlying_cluster)
        init_neighbors::Set{Int} = Set(
            collect(
                filter(neighbor -> !(neighbor in guarding_set), neighbors(underlying_cluster, vertex)),
            ),
        )
        vertices = [vertex]
        _grow_from_site(
            underlying_cluster,
            max_order,
            vertices,
            init_neighbors,
            guarding_set,
            out_array,
        )
        push!(guarding_set, vertex)
    end

    if !has_underlying_cluster(underlying_cluster)
        push!(out_array,
        Cluster(
            zeros(dimensions(underlying_cluster), 1),
            [1],
            zeros(Int, 2, 1, 1),
            zeros(1, 1),
            underlying_cluster,
            zeros(Int, 0),
            [zeros(Int, 0)],
            [zeros(Int, 0)],
            true,
            vertex_labeled(underlying_cluster),
            edge_labeled(underlying_cluster),
        ))
    end
    out_array
end

"""
Grows the subclusters from a specific site, up till specific order
and outputs them into the out_array. This adds all subclusters up
until the appropriate order, which is why it is only for the lattice
and not for individual clusters, which behave a bit differently

Inputs:
      lattice: cluster with coordinates as vertex labels,
      vertex colors, and edge weights

      max_order: Integer that details the order that
      the subclusters will be generated at

      subclusters_vertices: array of vertices that are the current
      subcluster

      neighbors: set of vertices that are neighbors to the current
      subcluster

      guarding_set: set of vertices not to visit

      out_array: array of subclusters of the underlying_cluster

Output:
      Technically, the output is has_int_leaf, but in practice, the
      output is the out_array that gets added to.
"""
function _grow_from_site(
    lattice::Cluster{false, <:Any, <:Any},
    max_order::Integer,
    subcluster_vertices::AbstractVector{V},
    current_neighbors::Set{V},
    guarding_set::Set{V},
    out_array::AbstractVector{<:AbstractCluster},
) where {V<:Integer}

    push!(out_array, Cluster(lattice, subcluster_vertices))

    if length(subcluster_vertices) == max_order
        return true
    end

    has_int_leaf = false
    new_guarding_set = copy(guarding_set)

    while !isempty(current_neighbors)
        neighbor = pop!(current_neighbors)
        append!(subcluster_vertices, neighbor)

        new_neighbors = copy(current_neighbors)

        for vertex in neighbors(lattice, neighbor)
            if (
                !(vertex in subcluster_vertices) &
                !(vertex in new_guarding_set) &
                !(vertex in new_neighbors)
            )

                push!(new_neighbors, vertex)
            end
        end

        if _grow_from_site(
            lattice,
            max_order,
            subcluster_vertices,
            new_neighbors,
            new_guarding_set,
            out_array,
        )
            pop!(subcluster_vertices)
            has_int_leaf = true
        else
            pop!(subcluster_vertices)
            return (has_int_leaf)
        end
        push!(new_guarding_set, neighbor)
        if (nv_sv(lattice) - length(new_guarding_set)) < max_order
            return (has_int_leaf)
        end
    end
    return (has_int_leaf)
end

"""
Grows clusters from the given vertices of the underlying_cluster. Calls the code above, but
for each order until the max order in the cluster
"""
function grow(underlying_cluster::Cluster{true, <:Any, <:Any})
    out_array::Vector{AbstractCluster} = Vector()

    for max_order = 1:(nv_sv(underlying_cluster)-1)
        append!(out_array, grow(underlying_cluster, max_order))
    end

    if nv_sv(underlying_cluster) != nv(underlying_cluster)
    # add in the single site for cluster expansions
    append!(out_array,
    repeat([Cluster(
        zeros(dimensions(underlying_cluster), 1),
        [1],
        zeros(Int, 2, 1, 1),
        zeros(1, 1),
        underlying_cluster,
        zeros(Int, 0),
        [zeros(Int, 0)],
        [zeros(Int, 0)],
        true,
        vertex_labeled(underlying_cluster),
        edge_labeled(underlying_cluster),
    )], nv(underlying_cluster)) )
    end

    out_array
end

"""
Grows the subclusters from a specific site, up till specific order
and outputs them into the out_array.

Inputs:
      underlying_cluster: Graph with coordinates as vertex labels,
      vertex colors, and edge weights

      subclusters_vertices: array of vertices that are the current
      subcluster

      neighbors: set of vertices that are neighbors to the current
      subcluster

      guarding_set: set of vertices not to visit

      out_array: array of subclusters of the underlying_cluster

Output:
      Technically, the output is has_int_leaf, but in practice, the
      output is the out_array that gets added to.
"""
function _grow_from_site(
    underlying_cluster::Cluster{true, <:Any, <:Any},
    max_order,
    subcluster_vertices::AbstractVector{V},
    current_neighbors::AbstractSet{V},
    guarding_set::AbstractSet{V},
    out_array::AbstractVector{<:AbstractCluster},
) where {V<:Integer}

    if length(subcluster_vertices) == max_order
        push!(out_array, Cluster(underlying_cluster, subcluster_vertices))
        return true
    end

    has_int_leaf = false
    new_guarding_set = copy(guarding_set)

    while !isempty(current_neighbors)
        neighbor = pop!(current_neighbors)
        append!(subcluster_vertices, neighbor)

        new_neighbors = copy(current_neighbors)

        for vertex in neighbors(underlying_cluster, neighbor)
            if (
                !(vertex in subcluster_vertices) &
                !(vertex in new_guarding_set) &
                !(vertex in new_neighbors)
            )

                push!(new_neighbors, vertex)
            end
        end

        if _grow_from_site(
            underlying_cluster,
            max_order,
            subcluster_vertices,
            new_neighbors,
            new_guarding_set,
            out_array,
        )
            pop!(subcluster_vertices)
            has_int_leaf = true
        else
            pop!(subcluster_vertices)
            return (has_int_leaf)
        end
        push!(new_guarding_set, neighbor)
        if (nv_sv(underlying_cluster) - length(new_guarding_set)) < max_order
            return (has_int_leaf)
        end
    end
    return (has_int_leaf)
end

end

