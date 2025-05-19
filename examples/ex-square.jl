"""
Example generating the clusters for the square lattice site expansion
"""
module ex1

using NLCE

# Setting the basis, only one atom per basis here
basis = [[0, 0]]

# Setting the primitive vectors
primitive_vectors = [[1, 0], [0, 1]]

# Nearest neighbors are just distance 1 away
neighbors = [1]

max_order = 8

# Choosing site expansion since this is a lattice upon which we can
# perform the site expansion
nlce_clusters, bundle = site_expansion_NLCE(basis, primitive_vectors, neighbors, max_order)

# Writing all the files to the corresponding folder, creating the folder
# if it does not exist
filepath = "examples/outputs/ex-square"
mkpath(filepath)
filename = filepath * "/square_nn_$(max_order).json"

# Output in standard format, JSON file
write_to_file(nlce_clusters, bundle, filename)
end
