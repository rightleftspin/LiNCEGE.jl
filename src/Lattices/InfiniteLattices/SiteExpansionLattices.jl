struct SiteExpansionLattice <: AbstractInfiniteLattice
        max_order::UInt8
        unit_cell::UnitCell
        coordinates::Matrix{Int}
        adj_matrix::Matrix{Int}
        neighbor_list::Vector{ExpansionVertices{Int}}
end

function SiteExpansionLattice(max_order::Int, unit_cell::UnitCell)
        @assert max_order > 0 "max_order must be a positive integer"

        coordinates = generate_coordinates(max_order, basis_size(unit_cell), dimension(unit_cell))
        adj_matrix, neighbor_list = generate_lattice_data(coordinates, unit_cell)

        return SiteExpansionLattice(UInt8(max_order), unit_cell, coordinates, adj_matrix, neighbor_list)
end

centers(lattice::SiteExpansionLattice) = find_centers(lattice.coordinates)
max_order(lattice::SiteExpansionLattice) = lattice.max_order
n_unique_sites(lattice::SiteExpansionLattice) = basis_size(lattice.unit_cell)

neighbors(lattice::SiteExpansionLattice, vs::ExpansionVertices) =
        union(ExpansionVertices(), lattice.neighbor_list[vs])

get_coordinates(lattice::SiteExpansionLattice) = shift_unit_cell(lattice.unit_cell, lattice.coordinates)
get_labels(lattice::SiteExpansionLattice) = lattice.coordinates[end, :]
get_site_colors(lattice::SiteExpansionLattice) = lattice.unit_cell.site_colors[lattice.coordinates[end, :]]
bond_matrix(lattice::SiteExpansionLattice) = lattice.adj_matrix
