# NLCE.jl

Julia implementation of the Numerical Linked Cluster Expansion method for
generic infinite or finite lattices.

To use this repository while under development use the following commands:

This section only needs to be run the first time you clone the repository

julia --project=.

using Pkg

Pkg.instantiate()

exit()

From here on out, you can just run

julia --project=. examples/ex-1.jl

or any other file in the repository.

# Overview of the algorithm

## Infinite Lattices

Starts by taking

## Finite Lattices

- Generate Lattice
- Find all translationally invariant clusters on the lattice
- Partition set of all translationally invariant clusters
- Find NLCE weights
