"""
A julia package for computing LINked Cluster Expansions on a GEneral geometry (LINCEGE).

The general pipeline is:
1. Define a unit cell and lattice geometry.
2. Generate all unique clusters up to a desired order using a cluster set and hasher.
3. Build an `Expansion` from those clusters.
4. Call `summation!` to populate the NLCE weights.
5. Export results with `write_to_json`.

See the [online documentation](https://rightleftspin.github.io/Lincege.jl) for worked
examples on the square lattice, Kagome lattice, and Pyrochlore unit cell.
"""
module Lincege

using Base.Threads
using LinearAlgebra
using NautyGraphs
using JSON

# Some limited utility functions for the rest of the algorithm
include("util.jl")

# Basic not-implemented functionality, taken from Graphs.jl
include("NI.jl")

# All the necessary code for constructing an Expansion
include("Vertices/Vertices.jl")
include("UnitCells/UnitCells.jl")
include("Lattices/Lattices.jl")
include("Hashers/Hashers.jl")
include("Clusters/Clusters.jl")
include("Expansions/Expansions.jl")
include("Writers/Writers.jl")

export AbstractVertices, LatticeVertices, ExpansionVertices,
        Bond, UnitCell, ExpansionBond, ExpansionUnitCell, image_unit_cell,
        SiteExpansionLattice, StrongClusterExpansionLattice, WeakClusterExpansionLattice,
        TranslationClusterSet, IsomorphicClusterSet, SymmetricClusterSet,
        clusters_from_lattice!, clusters_from_clusters!,
        Expansion, summation!, write_to_json

# Extra Physics Related Code, generally slow and not needed for basic Cluster Expansion construction
#include("Physics/Physics.jl")
end
