```@meta
CurrentModule = Lincege
```

# LiNCEGE.jl

This is the documentation for [LiNCEGE.jl](https://github.com/rightleftspin/LiNCEGE.jl), a
julia package for computing LINked Cluster Expansions on a GEneral geometry (LiNCEGE).

The registered Julia package and module are named `Lincege`; load it with `using Lincege`.

The general pipeline that LiNCEGE follows is:
1. Define a unit cell and lattice geometry.
2. Generate all unique clusters up to a desired order using a cluster set and hasher.
3. Build an `Expansion` from those clusters.
4. Call `summation!` to populate the NLCE weights.
5. Export results with `write_to_json`.

Worked examples are shown in the side bar. A full computation from lattice information 
to computed observables for the square lattice Ising model is shown in 
[Square Lattice Ising Model](examples/square_lattice.md), which serves as a good starting point for new users.

For advanced users, see [API Reference](api.md) for the complete documentation.
