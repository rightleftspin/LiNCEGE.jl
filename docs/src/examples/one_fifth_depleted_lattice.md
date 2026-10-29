# One Fifth Depleted Lattice

This example computes the linked cluster expansion on a one fifth depleted lattice up to
order 3.

## Setup

```@setup one-fifth-depleted
using Lincege
if !@isdefined(output_dir)
    output_dir = joinpath(pkgdir(Lincege), "examples", "outputs")
end
mkpath(output_dir)
```

Define the one-fifth-depleted unit cell:

```@example one-fifth-depleted
using Lincege

one_fifth_basis = [[0.0, 0.0], [0.0, 1.0], [0.0, 2.0],[1.0,1.0]]
one_fifth_pvecs = [[2.0, 1.0], [-1.0, 2.0]]
one_fifth_bonds = [Bond(1,2,[0,0],1),Bond(1,3,[0,-1],1),Bond(2,4,[0,0],1), Bond(3,4,[0,1],1),Bond(2,3,[0,0],2),Bond(1,4,[-1,0],2)]
one_fifth_uc = UnitCell(one_fifth_basis, one_fifth_pvecs, one_fifth_bonds, [1, 1, 1, 1])
```

Four basis sites and two tilted primitive vectors define the depleted square
geometry. Zero-offset bonds connect sites within the cell; the other offsets
connect neighboring cells. Bond types 1 and 2 distinguish two interaction
families, while equal site colors treat all sites as the same species.

## Imaging the unit cell and lattice

Load Plots and ColorSchemes to activate LiNCEGE.jl's plotting capabilities.

```@example one-fifth-depleted
using Plots, ColorSchemes
image_unit_cell(one_fifth_uc)
```

```@example one-fifth-depleted
image_lattice(one_fifth_uc)
```

## Building the lattice and clusters

The following code generates clusters up to order 3, starting with translational
clusters and reducing to isomorphic clusters.

```@example one-fifth-depleted
m_order = 3
lattice = SiteExpansionLattice(m_order, one_fifth_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
reduce_clusters!(iso_clusters, trans_clusters)
```

## Computing the expansion

The final NLCE expansion weights are computed recursively.

```@example one-fifth-depleted
expansion = Expansion(iso_clusters, lattice)
summation!(expansion, m_order)
```

## Writing to JSON

Lastly, write the results to a JSON file for future computations.

```julia
write_to_json(expansion, lattice, joinpath(output_dir, "one_fifth_depleted_lattice_site_expansion.json"))
```
