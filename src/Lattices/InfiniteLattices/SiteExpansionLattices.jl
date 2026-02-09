struct SiteExpansionLattice <: AbstractLattice
    max_order::UInt8
    unit_cell::UnitCell
    n_unique_sites::Int
    coordinates::AbstractMatrix{Int}
    neighbor_list::AbstractVector{LatticeVertices}
end

centers(lattice::SiteExpansionLattice) = lattice.centers
max_order(lattice::SiteExpansionLattice) = lattice.max_order
n_unique_sites(lattice::SiteExpansionLattice) = lattice.n_unique_sites

neighbors(lattice::SiteExpansionLattice, vs::LatticeVertices) = union(LatticeVertices(), lattice.neighbor_list[vs])
