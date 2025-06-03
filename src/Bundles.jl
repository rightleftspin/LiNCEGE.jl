
abstract type AbstractBundle end

# TODO: Documentation
"""
"""
mutable struct SiteExpansionBundle <: AbstractBundle

        "Underlying lattice for this NLCE expansion"
        lattice::Cluster
        "Coordinates for the lattice, in sorted order"
        coordinates::AbstractMatrix{<:Real}
        "Sublattice Coordinates for the lattice, in sorted order"
        sublattice_coordinates::AbstractMatrix{<:Integer}
        "Start points on the lattice super vertices for the NLCE expansion"
        start::AbstractVector{<:Integer}
        "Maximum order for the NLCE expansion"
        max_order::Integer
        "Hashing function for lattice embedding constant count"
        hashing_fxn::Any
        "Cluster info for all clusters invariant under the hashing function of this bundle"
        cluster_info::Any

        function SiteExpansionBundle(
                lattice::Cluster,
                coordinates::AbstractMatrix{<:Real},
                sublattice_coordinates::AbstractMatrix{<:Integer},
                start::AbstractVector{<:Integer},
                max_order::Integer,
                hashing_fxn,
        )

                @assert nv(lattice) == nsv(lattice) "Number of Vertices and Super Vertices need to be equivalent for site expansion"
                @assert nv(lattice) == size(coordinates, 1) "Number of coordinates needs to be the number of sites in the lattice"
                @assert nv(lattice) == size(sublattice_coordinates, 1) "Number of sublattice coordinates needs to be the number of sites in the lattice"
                @assert size(coordinates, 2) == (size(sublattice_coordinates, 2) - 1) "Sublattice coordinate or real space coordiante dimension is incorrect"
                @assert issorted(eachrow(coordinates)) "Coordinates must be sorted"
                @assert maximum(start) <= nsv(lattice) "Some start points are not within the lattice"

                bundle = new()
                bundle.lattice = lattice
                bundle.coordinates = coordinates
                bundle.sublattice_coordinates = sublattice_coordinates
                bundle.start = start
                bundle.max_order = max_order
                bundle.hashing_fxn = hashing_fxn
                bundle
        end
end

"""
Generates SiteExpansionBundle with lattice of size 2 * max_order + 1
"""
function SiteExpansionBundle(
        basis::AbstractVector{<:AbstractVector{<:Real}},
        primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
        neighbors::AbstractVector{<:Real},
        max_order::Integer,
        hashing_fxn;
        basis_labels::AbstractVector{<:Integer}=repeat([1], length(basis)),
)
        primitive_lattice, unrotated_lattice =
                generate_primitive_lattice(primitive_vectors, max_order)

        unsorted_coords = add_basis_coords(basis, primitive_lattice)
        unsorted_sublattice_coords = add_basis_sublattice(basis, unrotated_lattice)

        sort_perm_coords = sortperm(eachrow(unsorted_coords))

        coords = unsorted_coords[sort_perm_coords, :]
        sublattice_coords = unsorted_sublattice_coords[sort_perm_coords, :]
        labels =
                repeat(basis_labels, fld(size(unsorted_coords, 1), length(basis)))[sort_perm_coords]
        translation_labels =
                repeat(1:length(basis), fld(size(unsorted_coords, 1), length(basis)))[sort_perm_coords]

        adj_list = adj_list_from_coords(coords, neighbors)
        adj_matrices = adj_matrices_strong(
                coords,
                neighbors,
                labels,
                translation_labels,
                1:size(unsorted_coords, 1),
        )

        start = findfirst.(isapprox.(basis), (eachrow(coords),))

        lattice = Cluster(
                adj_list,
                [[i] for i = 1:size(coords, 1)],
                adj_matrices,
                length(unique(basis_labels)) != 1,
                length(neighbors) != 1,
        )

        SiteExpansionBundle(lattice, coords, sublattice_coords, start, max_order, hashing_fxn)

end

mutable struct StrongClusterExpansionBundle <: AbstractBundle

        "Underlying lattice for this NLCE expansion"
        lattice::Cluster
        "Coordinates for the lattice, in sorted order"
        coordinates::AbstractMatrix{<:Real}
        "Start points on the lattice super vertices for the NLCE expansion"
        start::AbstractVector{<:Integer}
        "Maximum order for the NLCE expansion"
        max_order::Integer
        "Hashing function for lattice embedding constant count"
        hashing_fxn::Any
        "Cluster info for all clusters invariant under the hashing function of this bundle"
        cluster_info::Any

        function StrongClusterExpansionBundle(
                lattice::Cluster,
                coordinates::AbstractMatrix{<:Real},
                start::AbstractVector{<:Integer},
                max_order::Integer,
                hashing_fxn,
        )

                @assert nv(lattice) == size(coordinates, 1) "Number of coordinates needs to be the number of sites in the lattice"
                @assert issorted(eachrow(coordinates)) "Coordinates must be sorted"
                @assert maximum(start) <= nsv(lattice) "Some start points are not within the lattice"

                bundle = new()
                bundle.lattice = lattice
                bundle.coordinates = coordinates
                bundle.start = start
                bundle.max_order = max_order
                bundle.hashing_fxn = hashing_fxn
                bundle
        end
end

"""
"""
function StrongClusterExpansionBundle(
        expansion_basis::AbstractVector{<:AbstractVector{<:Real}},
        struct_per_basis::AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}},
        expansion_labels::AbstractVector{<:AbstractVector{<:Integer}},
        expansion_primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
        expansion_neighbors::AbstractVector{<:Real},
        neighbors::AbstractVector{<:Real},
        max_order::Integer,
        hashing_fxn;
        basis_labels::AbstractVector{<:AbstractVector{<:Integer}}=[
                repeat([1], length(st)) for st in struct_per_basis
        ],
)
        lattice, unrotated_lattice =
                generate_primitive_lattice(expansion_primitive_vectors, max_order)

        exp_coords = add_basis_coords(expansion_basis, lattice)
        exp_sublattice_coords = add_basis_sublattice(expansion_basis, unrotated_lattice)

        unsorted_coords,
        unsorted_connections,
        unsorted_rev_connections,
        unsorted_labels,
        unsorted_translation_labels = hashing_lattice_coords(
                exp_coords,
                exp_sublattice_coords,
                struct_per_basis,
                basis_labels,
                expansion_labels,
        )

        exp_adj_list = adj_list_from_coords(exp_coords, expansion_neighbors)
        unsorted_adj_matrices = adj_matrices_strong(
                unsorted_coords,
                neighbors,
                unsorted_labels,
                unsorted_translation_labels,
                collect(Iterators.flatten(unsorted_rev_connections)),
        )
        start_points = findfirst.(isapprox.(expansion_basis), (eachrow(exp_coords),))

        sort_perm_coords = sortperm(eachrow(unsorted_coords))

        coords = unsorted_coords[sort_perm_coords, :]
        labels = unsorted_labels[sort_perm_coords]
        translation_labels = unsorted_translation_labels[sort_perm_coords]
        rev_connections = unsorted_rev_connections[sort_perm_coords]
        adj_matrices = unsorted_adj_matrices[:, sort_perm_coords, sort_perm_coords]

        connections::Vector{Vector{Int}} = []
        for unsorted_connection in unsorted_connections
                temp_connection = []
                for site in unsorted_connection
                        append!(temp_connection, findfirst(==(site), sort_perm_coords))
                end
                push!(connections, temp_connection)
        end

        lattice = Cluster(
                exp_adj_list,
                connections,
                adj_matrices,
                (length(unique(Iterators.flatten(basis_labels))) > 1),
                length(neighbors) != 1,
        )

        StrongClusterExpansionBundle(lattice, coords, start_points, max_order, hashing_fxn)

end

mutable struct WeakClusterExpansionBundle <: AbstractBundle

        "Underlying lattice for this NLCE expansion"
        lattice::Cluster
        "Coordinates for the lattice, in sorted order"
        coordinates::AbstractMatrix{<:Real}
        "Start points on the lattice super vertices for the NLCE expansion"
        start::AbstractVector{<:Integer}
        "Maximum order for the NLCE expansion"
        max_order::Integer
        "Hashing function for lattice embedding constant count"
        hashing_fxn::Any
        "Cluster info for all clusters invariant under the hashing function of this bundle"
        cluster_info::Any

        function WeakClusterExpansionBundle(
                lattice::Cluster,
                coordinates::AbstractMatrix{<:Real},
                start::AbstractVector{<:Integer},
                max_order::Integer,
                hashing_fxn,
        )

                @assert nv(lattice) == size(coordinates, 1) "Number of coordinates needs to be the number of sites in the lattice"
                @assert issorted(eachrow(coordinates)) "Coordinates must be sorted"
                @assert maximum(start) <= nsv(lattice) "Some start points are not within the lattice"

                bundle = new()
                bundle.lattice = lattice
                bundle.coordinates = coordinates
                bundle.start = start
                bundle.max_order = max_order
                bundle.hashing_fxn = hashing_fxn
                bundle
        end
end

"""
"""
function WeakClusterExpansionBundle(
        expansion_basis::AbstractVector{<:AbstractVector{<:Real}},
        struct_per_basis::AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}},
        expansion_labels::AbstractVector{<:AbstractVector{<:Integer}},
        expansion_primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
        expansion_neighbors::AbstractVector{<:Real},
        neighbors::AbstractVector{<:Real},
        max_order::Integer,
        hashing_fxn;
        basis_labels::AbstractVector{<:AbstractVector{<:Integer}}=[
                repeat([1], length(st)) for st in struct_per_basis
        ],
)
        lattice, unrotated_lattice =
                generate_primitive_lattice(expansion_primitive_vectors, max_order)

        exp_coords = add_basis_coords(expansion_basis, lattice)
        exp_sublattice_coords = add_basis_sublattice(expansion_basis, unrotated_lattice)

        unsorted_coords,
        unsorted_connections,
        unsorted_rev_connections,
        unsorted_labels,
        unsorted_translation_labels = hashing_lattice_coords(
                exp_coords,
                exp_sublattice_coords,
                struct_per_basis,
                basis_labels,
                expansion_labels,
        )

        exp_adj_list = adj_list_from_coords(exp_coords, expansion_neighbors)
        unsorted_adj_matrices = adj_matrices_weak(
                unsorted_coords,
                neighbors,
                unsorted_labels,
                unsorted_translation_labels,
                unsorted_rev_connections,
        )
        start_points = findfirst.(isapprox.(expansion_basis), (eachrow(exp_coords),))

        sort_perm_coords = sortperm(eachrow(unsorted_coords))

        coords = unsorted_coords[sort_perm_coords, :]
        labels = unsorted_labels[sort_perm_coords]
        translation_labels = unsorted_translation_labels[sort_perm_coords]
        rev_connections = unsorted_rev_connections[sort_perm_coords]
        adj_matrices = unsorted_adj_matrices[:, sort_perm_coords, sort_perm_coords]

        connections::Vector{Vector{Int}} = []
        for unsorted_connection in unsorted_connections
                temp_connection = []
                for site in unsorted_connection
                        append!(temp_connection, findfirst(==(site), sort_perm_coords))
                end
                push!(connections, temp_connection)
        end

        lattice = Cluster(
                exp_adj_list,
                connections,
                adj_matrices,
                (length(unique(Iterators.flatten(basis_labels))) > 1),
                length(neighbors) != 1,
        )

        WeakClusterExpansionBundle(lattice, coords, start_points, max_order, hashing_fxn)

end
