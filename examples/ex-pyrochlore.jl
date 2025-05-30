"""
Example generating the clusters for the square lattice site expansion
"""
module extriangularcluster

using NLCE

# Expansion basis is the center of each tiling piece
expansion_basis = [[1/2, sqrt(3)/4]]
# The structure around each expansion basis site, with reference to the
# expansion basis as the center
struct_per_basis = [[[-1/2, -sqrt(3)/4], [1/2, -sqrt(3)/4], [0, sqrt(3)/4]]]
# Labels to give each site in the struct per basis
expansion_labels = [[1, 1, 1]]
# Primitive vectors to expand the basis in the expansion
expansion_primitive_vectors = [[1, 0], [1/2, sqrt(3)/2]]
# Nearest neighbors to connect the sites of the expansion together
expansion_neighbors = [1]
# Per site factor is usually the number of unique expansion labels
# this allows you to get the final weights of the expansion in terms
# of per site weight
per_site_factor = 1
# neighbors on the underlying lattice
neighbors = [1]
# Number of sites in the expansion to generate,
# in this case, the number of squares
max_order = 5

# Choosing site expansion since this is a lattice upon which we can
# perform the site expansion
nlce_clusters, bundle = weak_cluster_expansion_NLCE(
                                                    expansion_basis,
                                                    struct_per_basis,
                                                    expansion_labels,
                                                    expansion_primitive_vectors,
                                                    expansion_neighbors,
                                                    neighbors,
                                                    per_site_factor,
                                                    max_order)

# Writing all the files to the corresponding folder, creating the folder
# if it does not exist
filepath = "examples/outputs/ex-triangular-cluster"
mkpath(filepath)
filename = filepath * "/triangular_cluster1_nn_$(max_order).json"

# Output in standard format, JSON file
write_to_file(nlce_clusters, bundle, filename)
end
