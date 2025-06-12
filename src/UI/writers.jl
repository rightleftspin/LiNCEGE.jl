"""
These are writers that write a set of clusters to a file for future use. There are many styles
of writing to disk, but the standard will be JSON files
"""

function write_to_file(
    nlce_output::AbstractDict{<:Cluster,Vector{<:Real}},
    bundle::AbstractBundle,
    filename::AbstractString,
)
    clusters = []
    for (cluster, mults) in nlce_output
        cluster_dict = Dict()
        cluster_dict["NLCE Order"] = nsv(cluster)
        cluster_dict["Number of Sites"] = nv(cluster)
        cluster_dict["Site Labels"] = all_vertex_labels(cluster)
        # TODO: cluster_dict["coordinates"] = get_coordinates(bundle, cluster)
        cluster_dict["Weighted Bonds"] = weighted_edge_list(cluster)
        cluster_dict["Multiplicities"] = mults

        push!(clusters, cluster_dict)
    end

    open(filename, "w") do io
        JSON3.pretty(io, clusters)
    end
end

"""
Fortran format specifies for each cluster the number of bonds. Then it has each
bond listed below it in its own line, tab spaced. Between each cluster, there is
an empty line and at the end of the file is multiplicity for each cluster at
every order. For example, with a maximum order of 4 on a square lattice,
the file for the clusters of order 3 would look like this:

start of file >>>
2
1	2	1
1	3	1

0 0 6 -38
<<< end of file

Here there is only one cluster, it has 2 bonds. The bonds connect
vertices 1 and 2, and vertices 1 and 3 both with weight 1. The cluster
has overall multiplicity 0 for both order 1 and 2, multiplicity 6 for
order 3 and -38 for order 4.

A general example would look like this:

start of file >>>
number_of_bonds
vertex1 vertex2 bond_weight
vertex1 vertex3 bond_weight
...

number_of_bonds
vertex1 vertex2 bond_weight
vertex1 vertex3 bond_weight
...

multiplicity_1 multiplicity_2 ...
multiplicity_1 multiplicity_2 ...
...
<<< end of file
"""
function write_to_file_fortran(
    nlce_output::AbstractDict{<:Cluster,Vector{<:Real}},
    bundle::AbstractBundle,
    filename::AbstractString,
    max_order::Integer,
)

    nlce_files = [open(filename * "_$(i).txt", "w") for i = 1:max_order]
    sorted_clusters = sort(collect(keys(nlce_output)), by = nv)

    for cluster in sorted_clusters
        edges = weighted_edge_list(cluster)
        write(nlce_files[nv(cluster)], "$(length(edges))\n")
        for edge in edges
            write(nlce_files[nv(cluster)], "$(join(edge, '\t'))\n")
        end
        write(nlce_files[nv(cluster)], "\n")
    end

    for cluster in sorted_clusters
        write(nlce_files[nv(cluster)], "$(join(nlce_output[cluster], ' '))\n")
    end

    close.(nlce_files)
end
