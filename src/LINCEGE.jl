"""
A julia package for computing LINked Cluster Expansions on a GEneral geometry (LINCEGE).
"""
module LINCEGE

using Base.Threads
using LinearAlgebra
using NautyGraphs
using JSON3

# Basic not-implemented functionality, taken from Graphs.jl
include("NI.jl")

# All the nessecary code for constructing an Expansion
include("Vertices/Vertices.jl")
include("UnitCells/UnitCells.jl")
include("Lattices/Lattices.jl")
include("Hashers/Hashers.jl")
include("Clusters/Clusters.jl")
include("Expansions/Expansions.jl")

export AbstractVertices, LatticeVertices, ExpansionVertices,
        Bond, UnitCell, ExpansionBond, ExpansionUnitCell, image_unit_cell,
        SiteExpansionLattice, StrongClusterExpansionLattice,
        TranslationClusterSet, IsomorphicClusterSet, SymmetricClusterSet,
        clusters_from_lattice!, clusters_from_clusters!,
        Expansion, summation!, write_to_json

# Extra Physics Related Code, generally slow and not needed for basic Cluster Expansion construction
#include("Physics/Physics.jl")
end
