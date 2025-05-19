"""
Example generating the data for a square lattice and 
writing to a file meant to be read by fortran code
"""
module ex1

using NLCE
using JLD

kagome_lattice = Dict("Expansion Basis" => [[0, sqrt(3)/3], [0, -sqrt(3)/3]],
                              "Struct Per Basis" => [[[0, -sqrt(3)/3], [1/2, sqrt(3)/6], [-1/2, sqrt(3)/6]], [[0, sqrt(3)/3], [1/2, -sqrt(3)/6], [-1/2, -sqrt(3)/6]]],
                              "Expansion Labels" => [[1, 2, 3], [1, 3, 2]],
                              "Expansion Primitive Vectors" => [[1, sqrt(3)], [1, -sqrt(3)]],
                              "Expansion Neighbors" => [2 * sqrt(3)/3],
                              )

neighbors = [1]
max_order = 12

@time begin
kagome_nlce_bundle = NLCE.WeakClusterExpansionBundle(
    kagome_lattice["Expansion Basis"],
    kagome_lattice["Struct Per Basis"],
    kagome_lattice["Expansion Labels"],
    kagome_lattice["Expansion Primitive Vectors"],
    kagome_lattice["Expansion Neighbors"],
    neighbors,
    max_order,
    NLCE.isomorphic_pruning,
)
end

@time begin
kagome_lattice_cluster_info = NLCE.lattice_constants!(
    kagome_nlce_bundle,
    (length(unique(Iterators.flatten(kagome_lattice["Expansion Labels"])))),
    single_site=true,
)
end

@time begin
NLCE.subclusters!(kagome_nlce_bundle, true)
end

@time begin
final_weights = NLCE.final_clusters(kagome_nlce_bundle, true)
end

num_sites = []
bond_lists = []
multiplicities = []
for (cluster, mults) in final_weights
    push!(num_sites, NLCE.nv(cluster))
    push!(bond_lists, NLCE.weighted_edge_list(cluster))
    push!(multiplicities, mults)
end

save("./outputs/Kagome_Lattice_$(max_order).jld",
     "num_sites", num_sites,
     "bond_lists", bond_lists,
     "multiplicities", multiplicities
     )

## Set the basis, since there is only one atom, it is at [0, 0]
#basis = [[0, 0]]
#
## Choose the primitive vectors, there are two on a square lattice
#primitive_vec = [[1, 0], [0, 1]]
#
## Choosing nearest neighbors (within distance 1 from each other)
#neighborhood = [1]
#
## Setting the maximum order
#max_order = 10
#
## Generating all the clusters using this information
#nlce_clusters = simple_NLCE(basis, primitive_vec, neighborhood, max_order)
#
#cluster_mults = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#for (cluster, mults) in nlce_clusters
#    cluster_mults[NLCE.nv(cluster)] += mults[findfirst(!=(0), mults)]
#end
#
#println(cluster_mults)
#
## Writing all the files to the corresponding folder, creating the folder
## if it does not exist
#filepath = "examples/outputs/ex-1/square_nn"
#mkpath(filepath)
#filename = filepath * "/square_nn"
#
## Write all the files in the "fortran" format
#write_to_file(nlce_clusters, filename)
#
#
#"""
#Fortran format specifies for each cluster the number of bonds. Then it has each
#bond listed below it in its own line, tab spaced. Between each cluster, there is
#an empty line and at the end of the file is multiplicity for each cluster at
#every order. For example, with a maximum order of 4 on a square lattice,
#the file for the clusters of order 3 would look like this:
#
#start of file >>>
#2
#1	2	1
#1	3	1
#
#0 0 6 -38
#<<< end of file
#
#Here there is only one cluster, it has 2 bonds. The bonds connect
#vertices 1 and 2, and vertices 1 and 3 both with weight 1. The cluster
#has overall multiplicity 0 for both order 1 and 2, multiplicity 6 for
#order 3 and -38 for order 4.
#
#A general example would look like this:
#
#start of file >>>
#number_of_bonds
#vertex1 vertex2 bond_weight
#vertex1 vertex3 bond_weight
#...
#
#number_of_bonds
#vertex1 vertex2 bond_weight
#vertex1 vertex3 bond_weight
#...
#
#multiplicity_1 multiplicity_2 ...
#multiplicity_1 multiplicity_2 ...
#...
#<<< end of file
#"""
end
