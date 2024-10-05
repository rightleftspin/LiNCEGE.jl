module NLCE

using NautyGraphs
using StaticArrays

# Add the relevant structs
include("structs/NLCELattices.jl")
include("structs/NLCEClusters.jl")

# Add the relevant helper functions
include("helpers/pruning.jl")
include("helpers/util.jl")
include("helpers/wrappers.jl")

# Add the basic pipeline
include("pipeline/grow.jl")
include("pipeline/prune.jl")
include("pipeline/propagate.jl")
include("pipeline/combine.jl")

export 
    simple_NLCE,
    write_to_file,
    write_to_file_fortran

end
