module NLCE

using LinearAlgebra
using Distances
using Rotations
using NautyGraphs
using JSON3 

# Add the relevant structs
include("Clusters.jl")
include("Bundles.jl")

# Add various methods for the structs
include("bundle_methods.jl")
include("cluster_methods.jl")
include("grow.jl")
include("hashing.jl")

# Add the relevant helper functions
include("util.jl")
include("wrappers.jl")

# Add Ising Model Simulation
include("ising.jl")

export SiteExpansionBundle,
    StrongClusterExpansionBundle,
    WeakClusterExpansionBundle,
    site_expansion_NLCE,
    simple_NLCE,
    write_to_file,
    write_to_file_fortran

end
