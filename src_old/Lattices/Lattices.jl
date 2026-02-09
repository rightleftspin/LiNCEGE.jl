module Lattices

using LinearAlgebra
using StaticArrays

import LINCEGE:
    _NI,
    Vertices.AbstractVertices,
    Vertices.ExpansionVertices,
    Vertices.LatticeVertices

abstract type AbstractLattice end
abstract type AbstractSiteExpansionLattice <: AbstractLattice end
abstract type AbstractClusterExpansionLattice <: AbstractLattice end
abstract type AbstractRandomLattice <: AbstractLattice end

centers(lattice::AbstractLattice) = _NI("centers")
max_order(lattice::AbstractLattice) = _NI("max_order")
neighbors(lattice::AbstractLattice, vs::AbstractVertices) = _NI("neighbors")
translation_form(lattice::AbstractLattice, vs::AbstractVertices) = _NI("translation_form")
is_weighted(lattice::AbstractLattice) = _NI("is_weighted")
get_labels(lattice::AbstractLattice, vs::AbstractVertices) = _NI("get_labels")
get_isomorphic_matrix(lattice::AbstractLattice, vs::AbstractVertices) = _NI("get_isomorphic_matrix")
n_unique_sites(lattice::AbstractLattice) = _NI("n_unique_sites")

include("util/misc.jl")
include("util/coordinate_constructors.jl")
include("util/neighbor_functions.jl")
include("util/graph_methods.jl")
include("util/translation_labeling.jl")

include("SiteExpansionLattices.jl")
include("ClusterExpansionLattices.jl")

export AbstractLattice,
    SiteExpansionLattice,
    ClusterExpansionLattice,
    max_order,
    centers,
    neighbors
end
