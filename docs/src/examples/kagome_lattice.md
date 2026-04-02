# Kagome Lattice

This example computes the linked cluster expansion on a kagome lattice up to
order 3.

## Setup

Define the three-site kagome unit cell:

```@example kagome
using LINCEGE

kagome_basis = [[0.0, 0.0], [1.0, 0.0], [0.5, sqrt(3) / 2]]
kagome_pvecs = [[2.0, 0.0], [1.0, sqrt(3)]]
kagome_bonds = [Bond(1, 2, [0, 0], 1), Bond(2, 3, [0, 0], 1), Bond(3, 1, [0, 0], 1),
        Bond(1, 2, [-1, 0], 1), Bond(1, 3, [0, -1], 1), Bond(2, 3, [1, -1], 1)]
kagome_uc = UnitCell(kagome_basis, kagome_pvecs, kagome_bonds, [1, 1, 1])
```

## Building the lattice and clusters

```@example kagome
m_order = 3
lattice = SiteExpansionLattice(m_order, kagome_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
clusters_from_clusters!(iso_clusters, trans_clusters)
```

## Computing the expansion

```@example kagome
expansion = Expansion(iso_clusters, lattice, m_order)
summation!(expansion, m_order)
expansion.weights
```

## Writing to JSON

```julia
write_to_json(expansion, lattice, iso_clusters, "kagome_lattice.json")
```
