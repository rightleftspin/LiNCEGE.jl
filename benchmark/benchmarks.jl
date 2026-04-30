using Lincege
using BenchmarkTools

SUITE = BenchmarkGroup()

basis_sq = [[0.0, 0.0]]
pvecs_sq = [[1.0, 0.0], [0.0, 1.0]]
bonds_sq = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
uc_sq = UnitCell(basis_sq, pvecs_sq, bonds_sq, [1])
lattice_sq = SiteExpansionLattice(9, uc_sq)
trans_sq = TranslationClusterSet(lattice_sq)
clusters_from_lattice!(trans_sq, lattice_sq)
iso_sq = IsomorphicClusterSet(lattice_sq)
clusters_from_clusters!(iso_sq, trans_sq)

SUITE["square"] = BenchmarkGroup()
SUITE["square"]["clusters_from_lattice"] = @benchmarkable begin
        cs = TranslationClusterSet($lattice_sq)
        clusters_from_lattice!(cs, $lattice_sq)
end
SUITE["square"]["clusters_from_clusters"] = @benchmarkable begin
        iso = IsomorphicClusterSet($lattice_sq)
        clusters_from_clusters!(iso, $trans_sq)
end
SUITE["square"]["Expansion"] = @benchmarkable Expansion($iso_sq, $lattice_sq)
SUITE["square"]["summation"] = @benchmarkable begin
        e = Expansion($iso_sq, $lattice_sq)
        summation!(e, 9)
end

basis_kag = [[0.0, 0.0], [1.0, 0.0], [0.5, sqrt(3) / 2]]
pvecs_kag = [[2.0, 0.0], [1.0, sqrt(3)]]
bonds_kag = [Bond(1, 2, [0, 0], 1), Bond(2, 3, [0, 0], 1), Bond(3, 1, [0, 0], 1),
        Bond(1, 2, [-1, 0], 1), Bond(1, 3, [0, -1], 1), Bond(2, 3, [1, -1], 1)]
uc_kag = UnitCell(basis_kag, pvecs_kag, bonds_kag, [1, 1, 1])
lattice_kag = SiteExpansionLattice(4, uc_kag)
trans_kag = TranslationClusterSet(lattice_kag)
clusters_from_lattice!(trans_kag, lattice_kag)
iso_kag = IsomorphicClusterSet(lattice_kag)
clusters_from_clusters!(iso_kag, trans_kag)

SUITE["kagome"] = BenchmarkGroup()
SUITE["kagome"]["clusters_from_lattice"] = @benchmarkable begin
        cs = TranslationClusterSet($lattice_kag)
        clusters_from_lattice!(cs, $lattice_kag)
end
SUITE["kagome"]["summation"] = @benchmarkable begin
        e = Expansion($iso_kag, $lattice_kag)
        summation!(e, 4)
end

# Pyrochlore Unit Cell — StrongClusterExpansionLattice
pyro_basis = [[[1 / 2, 1 / 2, 1 / 2], [1 / 2, -1 / 2, -1 / 2], [-1 / 2, 1 / 2, -1 / 2], [-1 / 2, -1 / 2, 1 / 2]]]
pyro_pvecs = [[2.0, 2.0, 0.0], [2.0, 0.0, 2.0], [0.0, 2.0, 2.0]]
pyro_lbonds = [
        ExpansionBond([1, 1], [1, 2], [0, 0, 0], 1),
        ExpansionBond([1, 1], [1, 3], [0, 0, 0], 1),
        ExpansionBond([1, 1], [1, 4], [0, 0, 0], 1),
        ExpansionBond([1, 2], [1, 3], [0, 0, 0], 1),
        ExpansionBond([1, 2], [1, 4], [0, 0, 0], 1),
        ExpansionBond([1, 3], [1, 4], [0, 0, 0], 1),
        ExpansionBond([1, 1], [1, 2], [0, 0, 1], 1),
        ExpansionBond([1, 1], [1, 3], [0, 1, 0], 1),
        ExpansionBond([1, 1], [1, 4], [1, 0, 0], 1),
        ExpansionBond([1, 2], [1, 3], [0, 1, -1], 1),
        ExpansionBond([1, 2], [1, 4], [1, 0, -1], 1),
        ExpansionBond([1, 3], [1, 4], [1, -1, 0], 1),
]
pyro_ebonds = [
        Bond(1, 1, [1, 0, 0], 1), Bond(1, 1, [0, 1, 0], 1), Bond(1, 1, [0, 0, 1], 1),
        Bond(1, 1, [0, 1, -1], 1), Bond(1, 1, [1, 0, -1], 1), Bond(1, 1, [1, -1, 0], 1),
]
uc_pyro = ExpansionUnitCell(pyro_basis, pyro_pvecs, pyro_lbonds, pyro_ebonds, [[1, 1, 1, 1]])
lattice_pyro = StrongClusterExpansionLattice(3, uc_pyro)
trans_pyro = TranslationClusterSet(lattice_pyro)
clusters_from_lattice!(trans_pyro, lattice_pyro)
iso_pyro = IsomorphicClusterSet(lattice_pyro)
clusters_from_clusters!(iso_pyro, trans_pyro)

SUITE["pyrochlore"] = BenchmarkGroup()
SUITE["pyrochlore"]["clusters_from_lattice"] = @benchmarkable begin
        cs = TranslationClusterSet($lattice_pyro)
        clusters_from_lattice!(cs, $lattice_pyro)
end
SUITE["pyrochlore"]["clusters_from_clusters"] = @benchmarkable begin
        iso = IsomorphicClusterSet($lattice_pyro)
        clusters_from_clusters!(iso, $trans_pyro)
end
SUITE["pyrochlore"]["summation"] = @benchmarkable begin
        e = Expansion($iso_pyro, $lattice_pyro)
        summation!(e, 3)
end

# Square Cluster — WeakClusterExpansionLattice
sq_cluster_basis = [[[-1 / 2, -1 / 2], [-1 / 2, 1 / 2], [1 / 2, -1 / 2], [1 / 2, 1 / 2]]]
sq_cluster_pvecs = [[1.0, 1.0], [1.0, -1.0]]
sq_cluster_lbonds = [
        ExpansionBond([1, 1], [1, 2], [0, 0], 1),
        ExpansionBond([1, 1], [1, 3], [0, 0], 1),
        ExpansionBond([1, 2], [1, 4], [0, 0], 1),
        ExpansionBond([1, 3], [1, 4], [0, 0], 1),
]
sq_cluster_ebonds = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
uc_sq_cluster = ExpansionUnitCell(sq_cluster_basis, sq_cluster_pvecs, sq_cluster_lbonds, sq_cluster_ebonds, [[1, 1, 1, 1]])
lattice_sq_cluster = WeakClusterExpansionLattice(4, uc_sq_cluster)
trans_sq_cluster = TranslationClusterSet(lattice_sq_cluster)
clusters_from_lattice!(trans_sq_cluster, lattice_sq_cluster)
iso_sq_cluster = IsomorphicClusterSet(lattice_sq_cluster)
clusters_from_clusters!(iso_sq_cluster, trans_sq_cluster)

SUITE["square_cluster"] = BenchmarkGroup()
SUITE["square_cluster"]["clusters_from_lattice"] = @benchmarkable begin
        cs = TranslationClusterSet($lattice_sq_cluster)
        clusters_from_lattice!(cs, $lattice_sq_cluster)
end
SUITE["square_cluster"]["clusters_from_clusters"] = @benchmarkable begin
        iso = IsomorphicClusterSet($lattice_sq_cluster)
        clusters_from_clusters!(iso, $trans_sq_cluster)
end
SUITE["square_cluster"]["summation"] = @benchmarkable begin
        e = Expansion($iso_sq_cluster, $lattice_sq_cluster)
        summation!(e, 4)
end
