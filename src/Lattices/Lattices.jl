module Lattices

using LinearAlgebra

import LINCEGE:
        _NI,
        Vertices.AbstractVertices,
        Vertices.ExpansionVertices,
        Vertices.LatticeVertices,
        UnitCells.UnitCell,
        UnitCells.dimension,
        UnitCells.basis_size,
        UnitCells.shift_unit_cell,
        UnitCells.find_possible_neighbors,
        UnitCells.neighbor_site

abstract type AbstractLattice end

centers(lattice::AbstractLattice) = _NI("centers")
max_order(lattice::AbstractLattice) = _NI("max_order")
neighbors(lattice::AbstractLattice, vs::AbstractVertices) = _NI("neighbors")
get_coordinates(lattice::AbstractLattice) = _NI("get_coordinates")
get_labels(lattice::AbstractLattice) = _NI("get_labels")
get_site_colors(lattice::AbstractLattice) = _NI("get_site_colors")
bond_matrix(lattice::AbstractLattice) = _NI("bond_matrix")
get_single_site_subgraphs(evs::ExpansionVertices, lattice::AbstractLattice) = _NI("get_single_site_subgraphs")

include("Connections.jl")
include("InfiniteLattices/InfiniteLattices.jl")

end
