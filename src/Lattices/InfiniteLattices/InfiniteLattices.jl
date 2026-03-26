abstract type AbstractInfiniteLattice <: AbstractLattice end

n_unique_sites(lattice::AbstractInfiniteLattice) = _NI("n_unique_sites")

include("util.jl")
include("SiteExpansionLattices.jl")
include("StrongClusterExpansionLattices.jl")
