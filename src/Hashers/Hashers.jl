module Hashers

using Base.Threads
using LinearAlgebra
using NautyGraphs

import LINCEGE:
        _NI,
        Vertices.ExpansionVertices,
        Lattices.centers,
        Lattices.max_order,
        Lattices.neighbors,
        Lattices.get_coordinates,
        Lattices.get_labels,
        Lattices.get_site_colors,
        Lattices.bond_matrix,
        Lattices.AbstractLattice,
        Lattices.SiteExpansionLattice

abstract type AbstractHasher end

ghash(h::AbstractHasher, evs::ExpansionVertices) = _NI("ghash")

include("util.jl")
include("TranslationHasher.jl")
include("IsomorphicHasher.jl")

end
