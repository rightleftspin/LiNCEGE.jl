"""
Example showing a cluster expansion on a square lattice
"""
module ex10

using NLCE

#sup_basis = [[0, 0, 0]]
#sub_basis = [[[0, 0, 0], [0, 1 / 4, 1 / 4], [1 / 4, 0, 1 / 4], [1 / 4, 1 / 4, 0]]]
#
#sup_primitive_vec = [[0, 1 / 2, 1 / 2], [1 / 2, 0, 1 / 2], [1 / 2, 1 / 2, 0]]
#sup_neighborhood = [sqrt(2) / 2]
#sub_neighborhood = [sqrt(2) / 4]

#sup_basis::Vector{Vector{Float64}} = [[0, 0]]
#sub_basis::Vector{Vector{Vector{Float64}}} = [[[0, 0], [1, 0], [0, 1], [1, 1]]]
#sup_primitive_vec::Vector{Vector{Float64}} = [[2, 0], [0, 2]]
#
#sup_neighborhood::Vector{Float64} = [2]
#sub_neighborhood::Vector{Float64} = [1]
#
#


max_order = 4
ssl_expansion_basis = [[0, 0], [2, 0]]
ssl_expansion_primitive_vectors = [[4, 0], [-2, 2]]
ssl_expansion_neighbors = [2]

ssl_struct_per_basis = [[[-0.5, 0], [0.5, 0]], [[0, 0.5], [0, -0.5]]]
ssl_colors = [[1, 2], [1, 2]]
ssl_neighbors = [1, sqrt(10) / 2]

println("start")

lattice = NLCE.Cluster(ssl_expansion_basis,
                       ssl_struct_per_basis,
                       ssl_expansion_primitive_vectors,
                       ssl_expansion_neighbors,
                       ssl_neighbors,
                       max_order)
println("lattice")

generated_clusters = NLCE.grow(lattice, max_order)
println("clusters")

iso_clusters = NLCE.prune(NLCE.isomorphic_pruning, Set(generated_clusters))
println("iso")

#for (hash, cluster_info) in iso_clusters
#    if NLCE.nv(cluster_info[1]) > 1
#        iso_clusters[hash] = (cluster_info[1], cluster_info[2] // 4, cluster_info[3])
#    end
#end

propogated = NLCE.propogate(NLCE.isomorphic_pruning, iso_clusters)
println("propogated")
for (hash, cluster) in propogated
    println("$(hash): $(cluster[1])")
    for (hash, mult) in cluster[3]
        println(hash)
        println(mult)
        end
end
println("--------------------------")


#sums = NLCE.nlce_summation(propogated, max_order)


# Initialize an empty output dictionary
output_dict = Dict{NLCE.Cluster, Vector{<:Real}}()

# Return the final sum for all clusters
for order = 1:(max_order + 1)
    for (hash, mult) in NLCE.nlce_summation(propogated, order)
        output_dict[iso_clusters[hash][1]] =
            append!(get(output_dict, iso_clusters[hash][1], Vector{Real}()), mult)
    end
end

# Writing all the files to the corresponding folder, creating the folder
# if it does not exist
filepath = "examples/outputs/ex-10/ssl_dimer"
mkpath(filepath)
filename = filepath * "/ssl_dimer"

# Write all the files in the default format
write_to_file(output_dict, filename)


# TODO: Deal with single site multiplicity in here this could be put anywhere tbh


## Generating all the clusters using this information
#nlce_clusters = simple_NLCE(basis, primitive_vec, neighborhood, max_order)
#
## Writing all the files to the corresponding folder, creating the folder
## if it does not exist
#filepath = "examples/outputs/ex-2/triangle_nn"
#mkpath(filepath)
#filename = filepath * "/triangle_nn"
#
## Write all the files in the "fortran" format
#write_to_file_fortran(nlce_clusters, filename, max_order)

end
