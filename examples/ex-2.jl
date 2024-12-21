"""
Example generating the data for a triangular lattice and 
writing to a file meant to be read by fortran code
"""
module ex2

using NLCE

# Set the basis, since there is only one atom, it is at [0, 0]
# Remember to list all these as floats
basis::Vector{Vector{Float64}} = [[0, 0]]

# Choose the primitive vectors, there are two on a square lattice
primitive_vec::Vector{Vector{Float64}} = [[0.5, sqrt(3) / 2], [1, 0]]

# Choosing nearest neighbors (within distance 1 from each other)
neighborhood::Vector{Float64} = [1]

# Setting the maximum order
max_order = 6

# Generating all the clusters using this information
nlce_clusters = simple_NLCE(basis, primitive_vec, neighborhood, max_order)

# Writing all the files to the corresponding folder, creating the folder
# if it does not exist
filepath = "examples/outputs/ex-2/triangle_nn"
mkpath(filepath)
filename = filepath * "/triangle_nn"

# Write all the files in the "fortran" format
write_to_file_fortran(nlce_clusters, filename, max_order)

end
