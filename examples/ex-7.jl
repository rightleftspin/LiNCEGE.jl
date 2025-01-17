"""
Example generating the data for a square lattice using symmetric
hashing alongside isomorphic hashing. Writing to a file with coordinates
of each site in a cluster
"""
module ex7

using NLCE

basis = [[0, 0]]

# Choose the primitive vectors, there are two on a square lattice
primitive_vec = [[1, 0], [0, 1]]

# Choosing nearest neighbors (within distance 1 from each other)
neighborhood = [1]

max_order = 4

# Generating all the clusters using this information
nlce_clusters = coord_NLCE(NLCE.square_symmetries, basis, primitive_vec, neighborhood, max_order)

# Writing all the files to the corresponding folder, creating the folder
# if it does not exist
filepath = "examples/outputs/ex-7/square_nn_sym"
mkpath(filepath)
filename = filepath * "/square_nn_sym"

# Write all the files in the default format
NLCE.write_to_file_coordinates(nlce_clusters, filename)

end
