using Profile

import LINCEGE:
    Lattices.SiteExpansionLattice,
    ClusterCollections.TranslationClusters,
    ClusterCollections.IsomorphicClusters,
    ClusterExpansions.SiteExpansion,
    ClusterExpansions.summation!

basis = [[0.0, 0.0]]
primitive_vectors = [[1, 0], [0, 1]]
colors = [1]
neighbor_distances = [1.0]
max_order_square = 8

square_lattice = SiteExpansionLattice(
    max_order_square,
    basis,
    primitive_vectors,
    colors,
    neighbor_distances,
)
translation_clusters_square = TranslationClusters(square_lattice)
isomorphic_clusters_square = IsomorphicClusters(translation_clusters_square, square_lattice)
square_expansion = SiteExpansion(isomorphic_clusters_square, square_lattice)

@time summation!(square_expansion, isomorphic_clusters_square, square_lattice)
#println(square_expansion)
@time summation!(square_expansion, isomorphic_clusters_square, square_lattice)
@time summation!(square_expansion, isomorphic_clusters_square, square_lattice)

@profile summation!(square_expansion, isomorphic_clusters_square, square_lattice)
Profile.print()
