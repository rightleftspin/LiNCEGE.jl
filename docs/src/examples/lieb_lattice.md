# Lieb Lattice

This example computes the linked cluster expansion on a Lieb lattice up to
order 3.

## Setup

```@setup lieb
using Lincege
if !@isdefined(output_dir)
    output_dir = joinpath(pkgdir(Lincege), "examples", "outputs")
end
mkpath(output_dir)
```

Define the three-site Lieb unit cell:

```@example lieb
using Lincege

lieb_basis = [[0.0, 0.0], [0.0, 1.0], [1, 0]]
lieb_pvecs = [[2.0, 0.0], [0.0, 2.0]]
lieb_bonds = [Bond(1,3,[1,0],1), Bond(1,2,[0,1],2)]
lieb_uc = UnitCell(lieb_basis, lieb_pvecs, lieb_bonds, [1, 1, 1])
```

The three-site basis contains a corner site and two edge sites. The offsets
`[1, 0]` and `[0, 1]` connect to sites in neighboring cells. Bond types 1 and 2
distinguish the horizontal and vertical interactions, while equal site colors
treat all sites as the same species.

## Imaging the unit cell and lattice

Load Plots and ColorSchemes to activate LiNCEGE.jl's plotting capabilities.

```@example lieb
using Plots, ColorSchemes
image_unit_cell(lieb_uc)
```

```@example lieb
image_lattice(lieb_uc)
```

## Building the lattice and clusters

The following code generates clusters up to order 3, starting with translational
clusters and reducing to isomorphic clusters.

```@example lieb
m_order = 3
lattice = SiteExpansionLattice(m_order, lieb_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
reduce_clusters!(iso_clusters, trans_clusters)
```

## Computing the expansion

The final NLCE expansion weights are computed recursively.

```@example lieb
expansion = Expansion(iso_clusters, lattice)
summation!(expansion, m_order)
```

## Writing to JSON

Lastly, write the results to a JSON file for future computations.

```julia
write_to_json(expansion, lattice, joinpath(output_dir, "lieb_lattice_site_expansion.json"))
```
