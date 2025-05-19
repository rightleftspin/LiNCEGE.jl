"""
Example generating the data for a triangular lattice and 
writing to a file meant to be read by fortran code
"""
module ex2

using NLCE

# Set the basis, since there is only one atom, it is at [0, 0]
basis = [[0, 0]]

# Choose the primitive vectors, there are two on a square lattice
primitive_vectors = [[0.5, sqrt(3) / 2], [1, 0]]

# Choosing nearest neighbors (within distance 1 from each other)
neighbors = [1]

# Setting the maximum order
max_order = 6

# Perform NLCE for this, use site expansion
nlce_clusters, bundle = site_expansion_NLCE(basis, primitive_vectors, neighbors, max_order)

# Writing all the files to the corresponding folder, creating the folder
# if it does not exist
filepath = "examples/outputs/ex-triangular"
mkpath(filepath)
filename = filepath * "/triangular_nn_6.json"

# Write all the files in the "fortran" format
#write_to_file_fortran(nlce_clusters, bundle, filename, max_order)
write_to_file(nlce_clusters, bundle, filename)

end
