using Distances

struct StrongClusterExpansionLattice <: AbstractClusterExpansionLattice
    max_order::UInt8
    tiling_units::AbstractVector{<:AbstractMatrix{Float64}}
    primitive_vectors::AbstractMatrix{Float64}
    centers::ExpansionVertices
    coordinates::AbstractMatrix{Float64}
    sublattice_coordinates::AbstractMatrix{Int}
    labels::AbstractVector{Int}
    translation_labels::AbstractVector{Int}
    n_unique_sites::Int
    neighbor_list::AbstractVector{ExpansionVertices}
    isomorphic_matrix::AbstractMatrix{Int}
    weighted::Bool
    translation_matrix::AbstractMatrix{Int}
    connections::AbstractVector{LatticeVertices}
end

function StrongClusterExpansionLattice(
    max_order::Int,
    tiling_units::AbstractVector{<:AbstractMatrix{<:Real}},
    primitive_vectors::AbstractMatrix{<:Real},
    colors::AbstractVector{<:AbstractVector{Int}},
    neighbor_function::Function,
)

    @info "Starting lattice construction for order $(max_order)"
    # Constructing the expansion lattice
    (primitive_coordinates, cartesian_coordinates) =
        generate_primitive_coordinates(primitive_vectors, max_order)

    basis = vcat(find_center.(tiling_units)...)
    if size(basis, 1) == 1
        exp_neighbor_fxn = distance_neighbor_function([norm(primitive_vectors[:, 1])])
    else
        exp_neighbor_fxn = distance_neighbor_function(Float64[sort(unique(pairwise(euclidean, basis, dims=2)))[2]])
    end

    exp_coordinates = add_basis_coords(basis, primitive_coordinates)
    exp_sublattice_coordinates = add_basis_sublattice(basis, cartesian_coordinates)
    neighbor_list = adj_mat_to_adj_list_exp(exp_neighbor_fxn(exp_coordinates, exp_sublattice_coordinates))

    centers = find_coordinate_indices_exp(exp_coordinates, basis)

    # Constructing the real lattice
    coordinates, sublattice_coordinates, connections = add_real_space_coords(tiling_units, exp_coordinates, exp_sublattice_coordinates)

    labeled_tiling_units = label_tiling_units(tiling_units, primitive_vectors)
    translation_labels = label_coordinates(sublattice_coordinates, labeled_tiling_units)
    labels = label_coordinates(sublattice_coordinates, colors)

    isomorphic_matrix::Matrix{Int} = neighbor_function(coordinates, sublattice_coordinates)

    pw_dir = pairwise_direction(coordinates)
    translation_matrix::Matrix{Int}, unique_dirs = unique_direction_indices(pw_dir, isomorphic_matrix)

    @info "Lattice construction completed"

    StrongClusterExpansionLattice(
        UInt8(max_order),
        tiling_units,
        primitive_vectors,
        centers,
        coordinates,
        sublattice_coordinates,
        labels,
        translation_labels,
        length(unique(translation_labels)),
        neighbor_list,
        isomorphic_matrix,
        length(unique(isomorphic_matrix)) > 2,
        2 .^ translation_matrix,
        connections,
    )

end

function StrongClusterExpansionLattice(
    max_order::Int,
    tiling_units::AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}},
    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    colors::AbstractVector{<:AbstractVector{Int}},
    neighbor_function::Function,
)

    StrongClusterExpansionLattice(
        max_order,
        [to_matrix(basis) for basis in tiling_units],
        to_matrix(primitive_vectors),
        colors,
        neighbor_function,
    )
end

function StrongClusterExpansionLattice(
    max_order::Int,
    tiling_units::AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}},
    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    colors::AbstractVector{<:AbstractVector{Int}},
    neighbor_distances::AbstractVector{<:Real},
)

    StrongClusterExpansionLattice(
        max_order,
        tiling_units,
        primitive_vectors,
        colors,
        distance_neighbor_function(neighbor_distances),
    )
end

centers(lattice::AbstractClusterExpansionLattice) = lattice.centers
max_order(lattice::AbstractClusterExpansionLattice) = lattice.max_order
neighbors(lattice::AbstractClusterExpansionLattice, vs::ExpansionVertices) = union(ExpansionVertices(), lattice.neighbor_list[vs])
neighbors(lattice::AbstractClusterExpansionLattice, vs::LatticeVertices) = LatticeVertices()
is_weighted(lattice::AbstractClusterExpansionLattice) = lattice.weighted
n_unique_sites(lattice::AbstractClusterExpansionLattice) = lattice.n_unique_sites
connections(lattice::AbstractClusterExpansionLattice, vs::ExpansionVertices) = lattice.connections[vs]

function translation_form(lattice::StrongClusterExpansionLattice, evs::ExpansionVertices)
    lvs = union(LatticeVertices(), lattice.connections[evs])
    sum(lattice.translation_matrix[lvs, lvs], dims=2) + lattice.translation_labels[lvs]
end

function translation_form(lattice::StrongClusterExpansionLattice, lvs::LatticeVertices)
    lattice.translation_labels[lvs]
end

function get_isomorphic_matrix(lattice::StrongClusterExpansionLattice, evs::ExpansionVertices)
    lvs = union(LatticeVertices(), lattice.connections[evs])
    lattice.isomorphic_matrix[lvs, lvs]
end

function get_labels(lattice::StrongClusterExpansionLattice, evs::ExpansionVertices)
    lvs = union(LatticeVertices(), lattice.connections[evs])
    lattice.labels[lvs]
end

function get_isomorphic_matrix(lattice::StrongClusterExpansionLattice, lvs::LatticeVertices)
    lattice.isomorphic_matrix[lvs, lvs]
end

function get_labels(lattice::StrongClusterExpansionLattice, lvs::LatticeVertices)
    lattice.labels[lvs]
end
