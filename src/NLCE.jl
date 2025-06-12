module NLCE

using LinearAlgebra
using Distances
using DataStructures
using StaticArrays
using NautyGraphs
using JSON3

# Add Adjacency Matrix Graphs
include("AdjMatrixGraphs/AdjMatrixGraphs.jl")
include("AdjMatrixGraphs/BondGraphs.jl")
include("AdjMatrixGraphs/DirectionGraphs.jl")
include("AdjMatrixGraphs/DistanceGraphs.jl")
include("AdjMatrixGraphs/ExpansionLatticeGraphs.jl")

# Add Clusters
include("Clusters/Clusters.jl")
include("Clusters/BondClusters.jl")
include("Clusters/DirectionClusters.jl")
include("Clusters/DistanceClusters.jl")

# Add Lattices
include("Lattices/RealSpaceLattices.jl")
include("Lattices/ExpansionLattices.jl")

# Include Cluster Tree Traversals
include("Traversal/dfs.jl")
include("Traversal/vsimple.jl")

#

export SiteExpansionBundle,
    StrongClusterExpansionBundle,
    WeakClusterExpansionBundle,
    site_expansion_NLCE,
    simple_NLCE,
    write_to_file,
    write_to_file_fortran

end
