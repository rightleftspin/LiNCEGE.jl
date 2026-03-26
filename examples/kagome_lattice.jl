import LINCEGE:
        UnitCells.UnitCell,
        UnitCells.Bond,
        Lattices.SiteExpansionLattice,
        Clusters.TranslationClusterSet,
        Clusters.IsomorphicClusterSet,
        Clusters.clusters_from_lattice!,
        Clusters.clusters_from_clusters!,
        Expansions.SiteExpansion,
        Expansions.summation!,
        Expansions.write_to_json

kagome_basis = [[0.0, 0.0], [1.0, 0.0], [0.5, sqrt(3) / 2]]
kagome_pvecs = [[2.0, 0.0], [1.0, sqrt(3)]]
kagome_bonds = [Bond(1, 2, [0, 0], 1), Bond(2, 3, [0, 0], 1), Bond(3, 1, [0, 0], 1),
        Bond(1, 2, [-1, 0], 1), Bond(1, 3, [0, -1], 1), Bond(2, 3, [1, -1], 1)]
kagome_uc = UnitCell(kagome_basis, kagome_pvecs, kagome_bonds, [1, 1, 1])

m_order = 3
lattice = SiteExpansionLattice(m_order, kagome_uc)
trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)
iso_clusters = IsomorphicClusterSet(lattice)
clusters_from_clusters!(iso_clusters, trans_clusters)

expansion = SiteExpansion(iso_clusters, lattice, m_order)
summation!(expansion, m_order)


path = "./examples/output/kagome_lattice.json"
write_to_json(expansion, lattice, iso_clusters, path)

