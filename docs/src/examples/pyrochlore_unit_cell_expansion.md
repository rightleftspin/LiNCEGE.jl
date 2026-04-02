# Pyrochlore Lattice Unit Cell Expansion

This example computes the linked cluster expansion on the pyrochlore lattice up
to order 3, this utilizes the unit cell expansion on the lattice.

## Setup

Define the unit cell with basis positions, primitive vectors, and
nearest-neighbour bonds:

```@example pyrochlore
using LINCEGE

pyro_exp_uc_basis = [
        [
                [1 / 2, 1 / 2, 1 / 2],
                [1 / 2, -1 / 2, -1 / 2],
                [-1 / 2, 1 / 2, -1 / 2],
                [-1 / 2, -1 / 2, 1 / 2],
        ],
]
pyro_exp_uc_pvecs = [[2.0, 2.0, 0.0], [2.0, 0.0, 2.0], [0.0, 2.0, 2.0]]
pyro_exp_uc_lbonds = [
        # Inside a Unit Cell
        ExpansionBond([1, 1], [1, 2], [0, 0, 0], 1),
        ExpansionBond([1, 1], [1, 3], [0, 0, 0], 1),
        ExpansionBond([1, 1], [1, 4], [0, 0, 0], 1),
        ExpansionBond([1, 2], [1, 3], [0, 0, 0], 1),
        ExpansionBond([1, 2], [1, 4], [0, 0, 0], 1),
        ExpansionBond([1, 3], [1, 4], [0, 0, 0], 1),
        # Between Unit Cells
        ExpansionBond([1, 1], [1, 2], [0, 0, 1], 1),
        ExpansionBond([1, 1], [1, 3], [0, 1, 0], 1),
        ExpansionBond([1, 1], [1, 4], [1, 0, 0], 1),
        ExpansionBond([1, 2], [1, 3], [0, 1, -1], 1),
        ExpansionBond([1, 2], [1, 4], [1, 0, -1], 1),
        ExpansionBond([1, 3], [1, 4], [1, -1, 0], 1),
]
pyro_exp_uc_ebonds = [
        Bond(1, 1, [1, 0, 0], 1),
        Bond(1, 1, [0, 1, 0], 1),
        Bond(1, 1, [0, 0, 1], 1),
        Bond(1, 1, [0, 1, -1], 1),
        Bond(1, 1, [1, 0, -1], 1),
        Bond(1, 1, [1, -1, 0], 1),
]
pyro_exp_uc_uc = ExpansionUnitCell(
        pyro_exp_uc_basis,
        pyro_exp_uc_pvecs,
        pyro_exp_uc_lbonds,
        pyro_exp_uc_ebonds,
        [[1, 1, 1, 1]]
)
```

## Building the lattice and clusters

Note that we use a StrongClusterExpansionLattice because the expansion does not
share sites between different expansion units. If we were to do the tetrahedra
expansion, we would use a WeakClusterExpansionLattice instead.

```@example pyrochlore
m_order = 3
lattice = StrongClusterExpansionLattice(m_order, pyro_exp_uc_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
clusters_from_clusters!(iso_clusters, trans_clusters)
```

## Computing the expansion

```@example pyrochlore
expansion = Expansion(iso_clusters, lattice, m_order)
summation!(expansion, m_order)
expansion.weights
```

## Writing to JSON

```julia
write_to_json(expansion, lattice, iso_clusters, "pyrochlore_lattice_uc_exp.json")
```
