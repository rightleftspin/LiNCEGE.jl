# Shastry-Sutherland Lattice

This example computes the linked cluster expansion on a Shastry-Sutherland lattice up to
order 3.

## Setup

```@setup shastry-sutherland
using Lincege
if !@isdefined(output_dir)
    output_dir = joinpath(pkgdir(Lincege), "examples", "outputs")
end
mkpath(output_dir)
```

Define the unit cell of the Shastry-Sutherland lattice:

```@example shastry-sutherland
using Lincege

shsu_basis = [[0.0, 0.0], [(sqrt(3) + 1) / (2 * sqrt(2)), (sqrt(3) - 1) / (2 * sqrt(2))], [sqrt(2) / 2, sqrt(6) / 2], [(sqrt(3) + 3) / (2 * sqrt(2)), (sqrt(3) + 1) / (2 * sqrt(2))]]
shsu_pvecs = [[(sqrt(3) + 1) / sqrt(2), 0.0], [0.0, (sqrt(3) + 1) / sqrt(2)]]
shsu_bonds = [Bond(1, 2, [0, 0], 1), Bond(2, 3, [0, 0], 1), Bond(1, 4, [-1, 0], 1), Bond(3, 4, [-1, 0], 1),
        Bond(3, 4, [0, 0], 2), Bond(2, 3, [0, -1], 2), Bond(1, 2, [-1, 0], 2), Bond(1, 4, [-1, -1], 2),
        Bond(2, 4, [-1, -1], 3), Bond(1, 3, [-1, 0], 3)]
shsu_uc = UnitCell(shsu_basis, shsu_pvecs, shsu_bonds, [1, 2, 3, 4])
```

The four-site basis uses three bond types to distinguish the interaction
families. Nonzero bond offsets connect to neighboring cells. The four distinct
site colors keep the basis sites distinguishable when comparing clusters.

## Imaging the unit cell and lattice

Load Plots and ColorSchemes to activate LiNCEGE.jl's plotting capabilities.

```@example shastry-sutherland
using Plots, ColorSchemes
image_unit_cell(shsu_uc)
```

```@example shastry-sutherland
image_lattice(shsu_uc)
```

## Building the lattice and clusters

The following code generates clusters up to order 3, starting with translational
clusters and reducing to isomorphic clusters.

```@example shastry-sutherland
m_order = 3
lattice = SiteExpansionLattice(m_order, shsu_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
reduce_clusters!(iso_clusters, trans_clusters)
```

## Computing the expansion

The final NLCE expansion weights are computed recursively.

```@example shastry-sutherland
expansion = Expansion(iso_clusters, lattice)
summation!(expansion, m_order)
```

## Writing to JSON

Lastly, write the results to a JSON file for future computations.

```julia
write_to_json(expansion, lattice, joinpath(output_dir, "shastry_sutherland_lattice_site_expansion.json"))
```
