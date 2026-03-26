struct StrongClusterExpansionLattice <: AbstractInfiniteLattice
        max_order::UInt8
        unit_cell::UnitCell
        lattice_unit_cells::Vector{UnitCell}
        coordinates::Matrix{Int}
        adj_matrix::Matrix{Int}
        neighbor_list::Vector{ExpansionVertices{Int}}
        connections::StrongClusterConnections
end

function StrongClusterExpansionLattice(max_order::Int, expansion_unit_cell::UnitCell, lattice_unit_cells::Vector{UnitCell})
        @assert max_order > 0 "max_order must be a positive integer"
        @assert length(lattice_unit_cells) == basis_size(expansion_unit_cell) "Need one lattice_unit_cell per expansion site type"

        return StrongClusterExpansionLattice(
                UInt8(max_order),
                expansion_unit_cell,
                lattice_unit_cells,
                phys_coords,
                phys_adj_matrix,
                expansion_neighbor_list,
                StrongClusterConnections(connection_list)
        )
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
