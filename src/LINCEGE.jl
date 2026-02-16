"""
A julia package for computing LINked Cluster Expansions on a GEneral geometry (LINCEGE).
"""
module LINCEGE

using Base.Threads
using LinearAlgebra
using Distances
using StaticArrays
using NautyGraphs
using JSON3

# Basic not-implemented functionality, taken from Graphs.jl
include("NI.jl")

# All the nessecary code for constructing an Expansion
include("Vertices/Vertices.jl")
include("UnitCells/UnitCells.jl")
include("Lattices/Lattices.jl")
include("Clusters/Clusters.jl")
#include("Expansions/Expansions.jl")

# Extra Physics Related Code, generally slow and not needed for basic Cluster Expansion construction
#include("Physics/Physics.jl")
end
