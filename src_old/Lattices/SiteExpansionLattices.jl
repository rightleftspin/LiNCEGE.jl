struct SiteExpansionLattice <: AbstractSiteExpansionLattice
    max_order::UInt8
    basis::AbstractMatrix{Float64}
    primitive_vectors::AbstractMatrix{Float64}
    centers::LatticeVertices
    coordinates::AbstractMatrix{Float64}
    sublattice_coordinates::AbstractMatrix{Int}
    labels::AbstractVector{Int}
    translation_labels::AbstractVector{Int}
    n_unique_sites::Int
    neighbor_list::AbstractVector{LatticeVertices}
    isomorphic_matrix::AbstractMatrix{Int}
    weighted::Bool
    translation_matrix::AbstractMatrix{Int}
end

function SiteExpansionLattice(
    max_order::Int,
    basis::AbstractMatrix{<:Real},
    primitive_vectors::AbstractMatrix{<:Real},
    colors::AbstractVector{Int},
    neighbor_function::Function,
)

    @info "Starting lattice construction for order $(max_order)"
    (primitive_coordinates, cartesian_coordinates) =
        generate_primitive_coordinates(primitive_vectors, max_order)

    coordinates = add_basis_coords(basis, primitive_coordinates)

    centers = find_coordinate_indices(coordinates, basis)

    sublattice_coordinates = add_basis_sublattice(basis, cartesian_coordinates)

    labels = Int[colors[i] for i in sublattice_coordinates[:, end]]

    isomorphic_matrix::Matrix{Int} = neighbor_function(coordinates, sublattice_coordinates)

    pw_dir = pairwise_direction(coordinates)
    translation_matrix::Matrix{Int}, unique_dirs = unique_direction_indices(pw_dir, isomorphic_matrix)

    neighbor_list::Vector{LatticeVertices} = adj_mat_to_adj_list(isomorphic_matrix)

    @info "Lattice construction completed"

    SiteExpansionLattice(
        UInt8(max_order),
        basis,
        primitive_vectors,
        centers,
        coordinates,
        sublattice_coordinates,
        labels,
        2 .^ (sublattice_coordinates[:, end] .+ length(unique_dirs)),
        size(basis, 1),
        neighbor_list,
        isomorphic_matrix,
        length(unique(isomorphic_matrix)) > 2,
        2 .^ translation_matrix,
    )

end

function SiteExpansionLattice(
    max_order::Int,
    basis::AbstractVector{<:AbstractVector{<:Real}},
    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    colors::AbstractVector{Int},
    neighbor_function::Function,
)

    SiteExpansionLattice(
        max_order,
        to_matrix(basis),
        to_matrix(primitive_vectors),
        colors,
        neighbor_function,
    )
end

function SiteExpansionLattice(
    max_order::Int,
    basis::AbstractVector{<:AbstractVector{<:Real}},
    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    colors::AbstractVector{Int},
    neighbor_distances::AbstractVector{<:Real},
)

    SiteExpansionLattice(
        max_order,
        basis,
        primitive_vectors,
        colors,
        distance_neighbor_function(neighbor_distances),
    )
end

centers(lattice::AbstractSiteExpansionLattice) = lattice.centers
max_order(lattice::AbstractSiteExpansionLattice) = lattice.max_order
neighbors(lattice::AbstractSiteExpansionLattice, vs::LatticeVertices) = union(LatticeVertices(), lattice.neighbor_list[vs])
translation_form(lattice::AbstractSiteExpansionLattice, vs::LatticeVertices) = sum(lattice.translation_matrix[vs, vs], dims=2) + lattice.translation_labels[vs]
is_weighted(lattice::AbstractSiteExpansionLattice) = lattice.weighted
get_labels(lattice::AbstractSiteExpansionLattice, vs::AbstractVertices) = lattice.labels[vs]
get_isomorphic_matrix(lattice::AbstractSiteExpansionLattice, vs::AbstractVertices) = lattice.isomorphic_matrix[vs, vs]
n_unique_sites(lattice::AbstractSiteExpansionLattice) = lattice.n_unique_sites
