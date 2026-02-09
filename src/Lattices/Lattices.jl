module Lattices

using LinearAlgebra
using StaticArrays

import LINCEGE:
    _NI,
    Vertices.AbstractVertices

abstract type AbstractLattice end

centers(lattice::AbstractLattice) = _NI("centers")
max_order(lattice::AbstractLattice) = _NI("max_order")
neighbors(lattice::AbstractLattice, vs::AbstractVertices) = _NI("neighbors")

end
