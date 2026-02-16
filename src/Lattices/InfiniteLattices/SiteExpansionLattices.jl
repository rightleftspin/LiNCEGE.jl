struct SiteExpansionLattice <: AbstractInfiniteLattice
    max_order::UInt8
    unit_cell::UnitCell
    coordinates::AbstractMatrix{<:Int}
    neighbor_list::AbstractVector{AbstractVertices}
end

function SiteExpansionLattice(max_order::Int, unit_cell::UnitCell)
    @assert max_order > 0 "max_order must be a positive integer"

    # Generate coordinates and neighbor list based on the unit cell and max order
    coordinates = generate_coordinates(max_order, basis_size(unit_cell), dimension(unit_cell))
    neighbor_list = generate_neighbor_list(coordinates, unit_cell)

    return SiteExpansionLattice(UInt8(max_order), unit_cell, coordinates, neighbor_list)
end

centers(lattice::SiteExpansionLattice) = find_centers(lattice.coordinates)
max_order(lattice::SiteExpansionLattice) = lattice.max_order
n_unique_sites(lattice::SiteExpansionLattice) = basis_size(lattice.unit_cell)

neighbors(lattice::SiteExpansionLattice, vs::ExpansionVertices) = union(ExpansionVertices(), lattice.neighbor_list[vs])

get_coordinates(lattice::SiteExpansionLattice) = shift_unit_cell(lattice.unit_cell, lattice.coordinates)
get_labels(lattice::SiteExpansionLattice) = lattice.coordinates[end, :]
bond_matrix(lattice::SiteExpansionLattice) = generate_weighted_adj_matrix(lattice.coordinates, lattice.unit_cell)
