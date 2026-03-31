using LINCEGE
using BenchmarkTools

SUITE = BenchmarkGroup()

basis_sq = [[0.0, 0.0]]
pvecs_sq = [[1.0, 0.0], [0.0, 1.0]]
bonds_sq = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
uc_sq = UnitCell(basis_sq, pvecs_sq, bonds_sq, [1])
lattice_sq = SiteExpansionLattice(6, uc_sq)
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
SUITE["square"]["Expansion"] = @benchmarkable Expansion($iso_sq, $lattice_sq, 6)
SUITE["square"]["summation"] = @benchmarkable begin
        e = Expansion($iso_sq, $lattice_sq, 6)
        summation!(e, 6)
end

basis_kag = [[0.0, 0.0], [1.0, 0.0], [0.5, sqrt(3) / 2]]
pvecs_kag = [[2.0, 0.0], [1.0, sqrt(3)]]
bonds_kag = [Bond(1, 2, [0, 0], 1), Bond(2, 3, [0, 0], 1), Bond(3, 1, [0, 0], 1),
        Bond(1, 2, [-1, 0], 1), Bond(1, 3, [0, -1], 1), Bond(2, 3, [0, 1], 1)]
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
        e = Expansion($iso_kag, $lattice_kag, 4)
        summation!(e, 4)
end
