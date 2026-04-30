# Square Lattice Cluster Expansion

This example computes the linked cluster expansion on the square lattice up to
order 3, this utilizes the corner-sharing square expansion

## Setup

Define the unit cell with basis positions, primitive vectors, and
nearest-neighbour bonds:

```@example square_cluster
using Lincege

square_cluster_basis = [[[-1 / 2, -1 / 2], [-1 / 2, 1 / 2], [1 / 2, -1 / 2], [1 / 2, 1 / 2]]]
square_cluster_pvecs = [[1.0, 1.0], [1.0, -1.0]]
square_cluster_lbonds = [
        ExpansionBond([1, 1], [1, 2], [0, 0], 1),
        ExpansionBond([1, 1], [1, 3], [0, 0], 1),
        ExpansionBond([1, 2], [1, 4], [0, 0], 1),
        ExpansionBond([1, 3], [1, 4], [0, 0], 1),
]
square_cluster_ebonds = [
        Bond(1, 1, [1, 0], 1),
        Bond(1, 1, [0, 1], 1),
]
square_cluster_uc = ExpansionUnitCell(
        square_cluster_basis,
        square_cluster_pvecs,
        square_cluster_lbonds,
        square_cluster_ebonds,
        [[1, 1, 1, 1]]
)
```

## Building the lattice and clusters

Note that we use a WeakClusterExpansionLattice because the expansion shares
sites between different expansion units (the corners of the squares).

```@example square_cluster
m_order = 3
lattice = WeakClusterExpansionLattice(m_order, square_cluster_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
clusters_from_clusters!(iso_clusters, trans_clusters)
```

## Computing the expansion

```@example square_cluster
expansion = Expansion(iso_clusters, lattice)
summation!(expansion, m_order)
```

## Writing to JSON

```julia
write_to_json(expansion, lattice, iso_clusters, "pyrochlore_lattice_uc_exp.json")
```
