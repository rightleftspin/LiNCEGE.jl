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

square_basis = [[0.0, 0.0]]
square_pvecs = [[1.0, 0.0], [0.0, 1.0]]
square_bonds = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
square_uc = UnitCell(square_basis, square_pvecs, square_bonds, [1])

m_order = 3
lattice = SiteExpansionLattice(m_order, square_uc)
trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)
iso_clusters = IsomorphicClusterSet(lattice)
clusters_from_clusters!(iso_clusters, trans_clusters)

expansion = SiteExpansion(iso_clusters, lattice, m_order)
summation!(expansion, m_order)

path = "./examples/output/square_lattice.json"
write_to_json(expansion, lattice, iso_clusters, path)
