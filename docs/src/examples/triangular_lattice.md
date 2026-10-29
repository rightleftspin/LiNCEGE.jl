# Triangular Lattice

This example computes the linked cluster expansion on a triangular lattice up to
order 3.

## Setup

```@setup triangular
using Lincege
if !@isdefined(output_dir)
    output_dir = joinpath(pkgdir(Lincege), "examples", "outputs")
end
mkpath(output_dir)
```

Define the triangular unit cell:

```@example triangular
using Lincege

triangular_basis = [[0.0, 0.0]]
triangular_pvecs = [[1.0, 0.0], [1.0/2, sqrt(3)/2]]
triangular_bonds = [Bond(1,1,[1,0],1),Bond(1,1,[0,1],1),Bond(1,1,[-1,1],1)]
triangular_uc = UnitCell(triangular_basis, triangular_pvecs, triangular_bonds, [1])
```

One site and two primitive vectors at 60 degrees generate the triangular
lattice. The offsets `[1, 0]`, `[0, 1]`, and `[-1, 1]` select three bond
directions; their reverses complete the six neighbors. Equal site colors
and bond types treat these sites and interactions as equivalent.

## Imaging the unit cell and lattice

Load Plots and ColorSchemes to activate LiNCEGE.jl's plotting capabilities.

```@example triangular
using Plots, ColorSchemes
image_unit_cell(triangular_uc)
```

```@example triangular
image_lattice(triangular_uc)
```

## Building the lattice and clusters

The following code generates clusters up to order 3, starting with translational
clusters and reducing to isomorphic clusters.

```@example triangular
m_order = 3
lattice = SiteExpansionLattice(m_order, triangular_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
reduce_clusters!(iso_clusters, trans_clusters)
```

## Computing the expansion

The final NLCE expansion weights are computed recursively.

```@example triangular
expansion = Expansion(iso_clusters, lattice)
summation!(expansion, m_order)
```

## Writing to JSON

Lastly, write the results to a JSON file for future computations.

```julia
write_to_json(expansion, lattice, joinpath(output_dir, "triangular_lattice_site_expansion.json"))
```
