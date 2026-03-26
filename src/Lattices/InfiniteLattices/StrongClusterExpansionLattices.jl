struct StrongClusterExpansionLattice <: AbstractInfiniteLattice
        max_order::UInt8
        unit_cell::UnitCell
        coordinates::Matrix{Int}
        adj_matrix::Matrix{Int}
        neighbor_list::Vector{ExpansionVertices{Int}}
        connections::Connections
end

function StrongClusterExpansionLattice(max_order::Int, unit_cell::UnitCell)
        @assert max_order > 0 "max_order must be a positive integer"

        coordinates = generate_coordinates(max_order, basis_size(unit_cell), dimension(unit_cell))
        adj_matrix, neighbor_list = generate_lattice_data(coordinates, unit_cell)

        return StrongClusterExpansionLattice(UInt8(max_order), unit_cell, coordinates, adj_matrix, neighbor_list)
end

centers(lattice::StrongClusterExpansionLattice) = find_centers(lattice.coordinates)
max_order(lattice::StrongClusterExpansionLattice) = lattice.max_order
n_unique_sites(lattice::StrongClusterExpansionLattice) = basis_size(lattice.unit_cell)

neighbors(lattice::StrongClusterExpansionLattice, vs::ExpansionVertices) =
        union(ExpansionVertices(), lattice.neighbor_list[vs])

get_coordinates(lattice::StrongClusterExpansionLattice) = shift_unit_cell(lattice.unit_cell, lattice.coordinates)
get_labels(lattice::StrongClusterExpansionLattice) = lattice.coordinates[end, :]
get_site_colors(lattice::StrongClusterExpansionLattice) = lattice.unit_cell.site_colors[lattice.coordinates[end, :]]
bond_matrix(lattice::StrongClusterExpansionLattice) = lattice.adj_matrix
connections(lattice::StrongClusterExpansionLattice) = lattice.connections
get_single_site_clusters(lattice::StrongClusterExpansionLattice) = []
get_single_site_subgraphs(lattice::StrongClusterExpansionLattice, evs::ExpansionVertices) = []
