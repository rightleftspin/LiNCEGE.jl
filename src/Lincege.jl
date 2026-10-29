"""
A julia package for computing LINked Cluster Expansions on a GEneral geometry (LiNCEGE).

The general pipeline is:
1. Define a unit cell and lattice geometry.
2. Generate all unique clusters up to a desired order using a cluster set and hasher.
3. Build an `Expansion` from those clusters.
4. Call `summation!` to populate the NLCE weights.
5. Export results with `write_to_json`.

See the [online documentation](https://rightleftspin.github.io/LiNCEGE.jl) for worked
examples on the square lattice, Kagome lattice, and Pyrochlore unit cell.
"""
module Lincege

using Base.Threads
using LinearAlgebra
using NautyGraphs
using JSON
using PrettyTables

# Basic not-implemented functionality, taken from Graphs.jl
include("NI.jl")

# All the necessary code for constructing an Expansion
include("Vertices/Vertices.jl")
include("UnitCells/UnitCells.jl")
# Utility helpers and example unit cells
include("util.jl")
include("Lattices/Lattices.jl")
include("Hashers/Hashers.jl")
include("Clusters/Clusters.jl")
include("Expansions/Expansions.jl")
include("Importers/Importers.jl")
include("Physics/Physics.jl")

# Vertices
export AbstractVertices, LatticeVertices, ExpansionVertices

# Unit Cells
export Bond, UnitCell, ExpansionBond, ExpansionUnitCell, image_unit_cell, image_lattice

# Lattices
export SiteExpansionLattice, StrongClusterExpansionLattice, WeakClusterExpansionLattice,
        FiniteLattice, FiniteStrongClusterExpansionLattice, FiniteWeakClusterExpansionLattice

# Hashers
export TranslationHasher, IsomorphicHasher, SymmetricHasher, ConnectedHasher

# Clusters
export TranslationClusterSet, IsomorphicClusterSet, SymmetricClusterSet, ConnectedClusterSet,
        clusters_from_lattice!, reduce_clusters!

# Expansions
export Expansion, summation!, weights, write_to_json, write_to_txt_file,
        print_latex_table, print_html_table, print_ascii_table,
        latex_table_column_labels

# Importers
export ImportedCluster, import_from_json, import_from_json_by_order

# Physics
export AbstractPhysicsSolver, observables, eigenvalues, eigenvectors, cluster_weights,
        IsingSolver,
        NLCEResult, perform_nlce,
        euler_resummation, wynn_resummation, apply_resummations
end
