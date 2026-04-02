```@meta
CurrentModule = LINCEGE
```

# LINCEGE

Documentation for [LINCEGE](https://github.com/rightleftspin/LINCEGE.jl), a
julia package for computing LINked Cluster Expansions on a GEneral geometry
(LINCEGE).

The general pipeline is:

    Define a unit cell and lattice geometry.
    Generate all unique clusters up to a desired order using a cluster set and hasher.
    Build an Expansion from those clusters.
    Call summation! to populate the NLCE weights.
    Export results with write_to_json.

See the online documentation for worked examples on the square lattice, Kagome
lattice, and Pyrochlore unit cell.

```@autodocs
Modules = [LINCEGE]
```
