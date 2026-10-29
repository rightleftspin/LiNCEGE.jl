# Kagome Lattice

This example computes the linked cluster expansion on a kagome lattice up to
order 3.

## Setup

```@setup kagome
using Lincege
if !@isdefined(output_dir)
    output_dir = joinpath(pkgdir(Lincege), "examples", "outputs")
end
mkpath(output_dir)
```

Define the three-site kagome unit cell:

```@example kagome
using Lincege

kagome_basis = [[0.0, 0.0], [1.0, 0.0], [0.5, sqrt(3) / 2]]
kagome_pvecs = [[2.0, 0.0], [1.0, sqrt(3)]]
kagome_bonds = [Bond(1, 2, [0, 0], 1), Bond(2, 3, [0, 0], 1), Bond(3, 1, [0, 0], 1),
        Bond(1, 2, [-1, 0], 1), Bond(1, 3, [0, -1], 1), Bond(2, 3, [1, -1], 1)]
kagome_uc = UnitCell(kagome_basis, kagome_pvecs, kagome_bonds, [1, 1, 1])
```

The kagome geometry consists of corner-sharing triangles with three bonds inside
the triangle and three outside. The offsets `[-1, 0]`, `[0, -1]`, and `[1, -1]` 
connect to neighboring triangles. Keeping all site colors and bond types equal
treats these sites and interactions as equivalent.

## Imaging the unit cell and lattice

Loading the Plots dependency automatically activates LiNCEGE.jl's plotting capabilities.

```@example kagome
using Plots, ColorSchemes
image_unit_cell(kagome_uc)
```

```@example kagome
image_lattice(kagome_uc)
```

## Building the lattice and clusters

The following code generates clusters up till order 3, starting with translational
clusters and reducing to isomorphic clusters.

```@example kagome
m_order = 3
lattice = SiteExpansionLattice(m_order, kagome_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
reduce_clusters!(iso_clusters, trans_clusters)
```

## Computing the expansion

The NLCE expansion final weights are computed recursively.

```@example kagome
expansion = Expansion(iso_clusters, lattice)
summation!(expansion, m_order)
```

## Writing to JSON

Lastly, the results are written to a JSON file to be used for future computations.

```julia
write_to_json(expansion, lattice, joinpath(output_dir, "kagome_lattice_site_expansion.json"))
```
