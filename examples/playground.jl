module playground

using NLCE

# Diamond Lattice
#basis = [[0, 0, 0], [1 / 4, 1 / 4, 1 / 4]]
#
#primitive_vec = [[0, 1 / 2, 1 / 2], [1 / 2, 0, 1 / 2], [1 / 2, 1 / 2, 0]]
#
#neighborhood = [sqrt(3) / 4]

# Pyrochlore Lattice
basis = [[0, 0, 0], [0, 1 / 4, 1 / 4], [1 / 4, 0, 1 / 4], [1 / 4, 1 / 4, 0]]

primitive_vec = [[0, 1 / 2, 1 / 2], [1 / 2, 0, 1 / 2], [1 / 2, 1 / 2, 0]]

neighborhood = [sqrt(2) / 4]

# Square Lattice
#basis = [[0, 0]]
#
#primitive_vec = [[1, 0], [0, 1]]
#
#neighborhood = [1]

# Triangular Lattice
#basis::Vector{Vector{Float64}} = [[0, 0]]
#
#primitive_vec::Vector{Vector{Float64}} = [[0.5, sqrt(3) / 2], [1, 0]]
#
#neighborhood::Vector{Float64} = [1]

max_order = 1

basis = basis - (sum(basis)./length(basis))

for sym in NLCE.pyrochlore_symmetries
    basis_new = []
    for elem in basis
        new_elem = sym * elem
        is_there = findfirst(isapprox(new_elem), basis)
        append!(basis_new, new_elem)
    end
end

file_name = "pyrochlore_nn_sym"

#nlce_clusters = coord_NLCE(NLCE.pyrochlore_symmetries, basis, primitive_vec, neighborhood, max_order)
#
## Writing all the files to the corresponding folder, creating the folder
## if it does not exist
#filepath = "examples/outputs/ex-5/" * file_name
#mkpath(filepath)
#filename = filepath * "/" * file_name

## Write all the files in the default format
#NLCE.write_to_file_coordinates(nlce_clusters, filename)

end
