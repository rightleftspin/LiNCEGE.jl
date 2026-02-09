module Lattices

abstract type AbstractInfiniteLattice end

num_unique_sites(lattice::AbstractInfiniteLattice) = _NI("n_unique_sites")

include("SiteExpansionLattice.jl")

export AbstractInfiniteLattice,
    num_unique_sites
end
