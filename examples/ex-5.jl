"""
Example generating the data for a pyrochlore lattice using symmetric
hashing alongside isomorphic hashing. Writing to a file with coordinates
"""
module ex5

using NLCE

# Diamond Lattice
#basis = [[0, 0, 0], [1 / 4, 1 / 4, 1 / 4]]
#
#primitive_vec = [[0, 1 / 2, 1 / 2], [1 / 2, 0, 1 / 2], [1 / 2, 1 / 2, 0]]
#
#neighborhood = [sqrt(3) / 4]

basis = [[0, 0, 0], [0, 1 / 4, 1 / 4], [1 / 4, 0, 1 / 4], [1 / 4, 1 / 4, 0]]

primitive_vec = [[0, 1 / 2, 1 / 2], [1 / 2, 0, 1 / 2], [1 / 2, 1 / 2, 0]]

neighborhood = [sqrt(2) / 4]

max_order = 3

# Generating all the clusters using this information
nlce_clusters, cluster_hashes, cluster_perms = coord_NLCE(NLCE.pyrochlore_symmetries, basis, primitive_vec, neighborhood, max_order)

# Writing all the files to the corresponding folder, creating the folder
# if it does not exist
filepath = "examples/outputs/ex-5/pyrochlore_nn_sym"
mkpath(filepath)
filename = filepath * "/pyrochlore_nn_sym_$(max_order)"

# Write all the files in the default format
NLCE.write_to_file_coordinates(nlce_clusters, cluster_hashes, cluster_perms, filename)

end
