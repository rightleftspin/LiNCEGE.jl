using Profile
using InteractiveUtils

import LINCEGE:
        UnitCells.UnitCell,
        UnitCells.Bond,
        Lattices.SiteExpansionLattice,
        Clusters.ClusterSet,
        Clusters.IsomorphicClusterSet,
        Clusters.clusters_from_lattice!,
        Clusters.clusters_from_clusters!,
        Expansions.SiteExpansion,
        Expansions.summation!
        #Expansions.faster_summation!,
        #Expansions.slow_summation!

basis = [[0.0, 0.0]]
primitive_vectors = [[1.0, 0.0], [0.0, 1.0]]
bonds = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
site_colors = [1]
unit_cell = UnitCell(basis, primitive_vectors, bonds, site_colors)

m_order = 12
lattice = SiteExpansionLattice(m_order, unit_cell)

translation_clusters = ClusterSet(lattice)
clusters_from_lattice!(translation_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
clusters_from_clusters!(iso_clusters, translation_clusters)

println("Summation Started")
expansion = SiteExpansion(iso_clusters, lattice, m_order)
@time summation!(expansion, m_order)
expansion = SiteExpansion(iso_clusters, lattice, m_order)
@time summation!(expansion, m_order)
expansion = SiteExpansion(iso_clusters, lattice, m_order)
@time summation!(expansion, m_order)
#println(expansion.weights)

#println("Slow Summation Started")
#expansion = SiteExpansion(iso_clusters, lattice, m_order)
#@time slow_summation!(expansion, m_order)
#expansion = SiteExpansion(iso_clusters, lattice, m_order)
#@time slow_summation!(expansion, m_order)
#expansion = SiteExpansion(iso_clusters, lattice, m_order)
#@time slow_summation!(expansion, m_order)
#println(expansion.weights)

#expansion = SiteExpansion(iso_clusters, lattice, m_order)
#@profile faster_summation!(expansion, m_order)
#Profile.print()
#
#expansion = SiteExpansion(iso_clusters, lattice, m_order)
#@profile slow_summation!(expansion, m_order)
#Profile.print()
#
#expansion = SiteExpansion(iso_clusters, lattice, m_order)
#@code_warntype slow_summation!(expansion, m_order)
