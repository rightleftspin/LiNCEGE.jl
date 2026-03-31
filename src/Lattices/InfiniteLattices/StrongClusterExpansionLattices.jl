struct StrongClusterExpansionLattice <: AbstractClusterExpansionLattice
        max_order::UInt8
        expansion_unit_cell::ExpansionUnitCell
        lattice_coordinates::Vector{Matrix{Int}}
        adj_matrix::Matrix{Int}
        neighbor_list::Vector{ExpansionVertices{Int}}
        connections::StrongClusterConnections
end

function StrongClusterExpansionLattice(max_order::Int, expansion_unit_cell::ExpansionUnitCell)
        @assert max_order > 0 "max_order must be a positive integer"

        return StrongClusterExpansionLattice(
                UInt8(max_order),
                expansion_unit_cell,
                expansion_coordinates,
                lattice_coordinates,
                adj_matrix,
                neighbor_list,
                StrongClusterConnections(connections_vec)
        )
end

centers(lattice::StrongClusterExpansionLattice) = ExpansionVertices(find_centers(lattice.expansion_coordinates))
max_order(lattice::StrongClusterExpansionLattice) = lattice.max_order
n_unique_sites(lattice::StrongClusterExpansionLattice) = sum(basis_size.(lattice.lattice_unit_cells))
n_site_colors(lattice::StrongClusterExpansionLattice) = length(unique(vcat([uc.site_colors for uc in lattice.lattice_unit_cells]...)))
neighbors(lattice::StrongClusterExpansionLattice, vs::ExpansionVertices) = union(ExpansionVertices(), lattice.neighbor_list[vs])
get_coordinates(lattice::StrongClusterExpansionLattice) = hcat([shift_unit_cell(uc, cs) for (uc, cs) in zip(lattice.lattice_unit_cells, lattice.lattice_coordinates)]...)
get_labels(lattice::StrongClusterExpansionLattice) = vcat([cs[end, :] for cs in lattice.lattice_coordinates]...)
get_site_colors(lattice::StrongClusterExpansionLattice) = vcat([uc.site_colors[cs[end, :]] for (uc, cs) in zip(lattice.lattice_unit_cells, lattice.lattice_coordinates)]...)
bond_matrix(lattice::StrongClusterExpansionLattice) = lattice.adj_matrix
connections(lattice::StrongClusterExpansionLattice) = lattice.connections
