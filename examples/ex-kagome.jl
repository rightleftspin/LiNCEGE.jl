"""
Example generating the data for a kagome lattice site expansion
along with the triangle expansion.
"""
module ex3

using NLCE

# All for the site expansion
basis = [[0, 0], [1, 0], [1 / 2, sqrt(3) / 2]]
primitive_vectors = [[2, 0], [1, sqrt(3)]]
neighbors = [1]
max_order = 5

nlce_clusters, bundle = site_expansion_NLCE(basis, primitive_vectors, neighbors, max_order)

# Writing all the files to the corresponding folder, creating the folder
# if it does not exist
filepath = "examples/outputs/ex-kagome"
mkpath(filepath)
filename = filepath * "/kagome_nn_$(max_order).json"

# Output in standard format, JSON file
# Be careful, this outputs multiplicities to floats since
# JSON cannot easily handle julia's rational number format
write_to_file(nlce_clusters, bundle, filename)

kagome_lattice_cluster_expansion = Dict(
        "Expansion Basis" => [[0, sqrt(3) / 3], [0, -sqrt(3) / 3]],
        "Struct Per Basis" => [
                [[0, -sqrt(3) / 3], [1 / 2, sqrt(3) / 6], [-1 / 2, sqrt(3) / 6]],
                [[0, sqrt(3) / 3], [1 / 2, -sqrt(3) / 6], [-1 / 2, -sqrt(3) / 6]],
        ],
        "Expansion Labels" => [[1, 2, 3], [1, 3, 2]],
        "Expansion Primitive Vectors" => [[1, sqrt(3)], [1, -sqrt(3)]],
        "Expansion Neighbors" => [2 * sqrt(3) / 3],
)
end
