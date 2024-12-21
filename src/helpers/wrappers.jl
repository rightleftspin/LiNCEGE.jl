"""
These are some high level convenience wrappers for the NLCE process.
"""

"""
Performs simple NLCE for most use cases.

Inputs: 
      basis: basis for the unit cell, this wrapper
      assumes all atoms in the basis are the same

      primitive_vectors: primitive vectors to tile the lattice

      neighborhood: array of distances to consider as neighbors. ie. The
      first nearest neighbor will be all neighbors with distance equal 
      (up to machine precision) to the first entry of this array

      max_order: highest NLCE order to go to

Output:
    hashmap containing hashes of clusters and their corresponding multiplicities
"""
function simple_NLCE(
    basis::AbstractVector{<:AbstractVector{<:Real}},
    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    neighborhood::AbstractVector{<:Real},
    max_order::Integer,
)

    # Create the lattice
    lattice = NLCELattice(basis, primitive_vectors, neighborhood, max_order)
    # Generate clusters on the lattice
    generated_clusters = grow(lattice, max_order)
    # find all the isomorphic clusters
    iso_clusters =
        prune(isomorphic_pruning, filtering(translational_pruning, generated_clusters))

    # Account for the size of the unit cell in the pruning
    for (hash, (cluster, mult, subcluster_mult)) in iso_clusters
        if nv(cluster) > 1
            iso_clusters[hash] = (cluster, mult // length(basis), subcluster_mult)
        end
    end

    # Find all their subclusters
    subclusters = propogate(isomorphic_pruning, iso_clusters)

    # Initialize an empty output dictionary
    output_dict = Dict{AbstractNLCECluster,Vector{<:Real}}()

    # Return the final sum for all clusters
    for order = 1:max_order
        for (hash, mult) in nlce_summation(subclusters, order)
            output_dict[iso_clusters[hash][1]] =
                append!(get(output_dict, iso_clusters[hash][1], Vector{Real}()), mult)
        end
    end

    output_dict
end

"""
Performs coordinate NLCE, this is useful if your observables depend on real distance.

Inputs: 
      basis: basis for the unit cell, this wrapper
      assumes all atoms in the basis are the same

      primitive_vectors: primitive vectors to tile the lattice

      neighborhood: array of distances to consider as neighbors. ie. The
      first nearest neighbor will be all neighbors with distance equal 
      (up to machine precision) to the first entry of this array

      max_order: highest NLCE order to go to

Output:
    hashmap containing clusters and their corresponding multiplicities
"""
function coord_NLCE(
    symmetries,
    basis::AbstractVector{<:AbstractVector{<:Real}},
    primitive_vectors::AbstractVector{<:AbstractVector{<:Real}},
    neighborhood::AbstractVector{<:Real},
    max_order::Integer,
)

    # Create the lattice
    lattice = NLCELattice(basis, primitive_vectors, neighborhood, max_order)
    # Generate clusters on the lattice
    generated_clusters = grow(lattice, max_order)
    # find all the symmetrically distinct clusters
    coordinates, colors, centers = generate_coordinates(
        basis,
        primitive_vectors,
        max_order,
        repeat([1], length(basis)),
    )

    permutations = find_permutations(coordinates, centers, symmetries)


    symmetric_pruning(cluster_pruning) = (
        hash(
            sort(
                NLCE.translational_form.([
                    cluster(
                        underlying_lattice(cluster_pruning),
                        perm[vertices(cluster_pruning)],
                    ) for perm in permutations
                ]),
            ),
        ),
        nothing,
    )

    sym_clusters =
        prune(symmetric_pruning, filtering(translational_pruning, generated_clusters))

    for (k, v) in sym_clusters
        println(k)
        println(v[2])
    end
    #iso_clusters = prune(isomorphic_pruning, sym_clusters)
    #count = 0
    #for (cluster, vals) in sym_clusters
    #    if nv(vals[1]) == 5
    #        count += vals[2]
    #    end
    #    println(nv(vals[1]))
    #end
    #println("")
    #println(count)

    # Account for the size of the unit cell in the pruning
    for (hash, (cluster, mult, subcluster_mult)) in sym_clusters
        if nv(cluster) > 1
            sym_clusters[hash] = (cluster, mult // length(basis), subcluster_mult)
        end
    end

    # Find all their subclusters
    subclusters = propogate(symmetric_pruning, sym_clusters)

    #for (hash, clusters) in subclusters


    # Initialize an empty output dictionary
    output_dict = Dict{AbstractNLCECluster,Vector{<:Real}}()

    # Return the final sum for all clusters
    for order = 1:max_order
        for (hash, mult) in nlce_summation(subclusters, order)
            output_dict[sym_clusters[hash][1]] =
                append!(get(output_dict, sym_clusters[hash][1], Vector{Real}()), mult)
        end
    end

    output_dict
end

function write_to_file_coordinates(
    nlce_output::AbstractDict{AbstractNLCECluster,Vector{<:Real}},
    filename::AbstractString,
)

    nlce_file = open(filename, "w")

    for (cluster, mults) in nlce_output
        write(nlce_file, "$(nv(cluster)):")
        for edge in edge_list(cluster)
            write(nlce_file, " $(join(edge, ' '))")
        end
        for coord in all_coordinates(cluster)
            write(nlce_file, " ($(join(coord, ',')))")
        end
        write(nlce_file, " : $(join(mults, ' '))\n")
    end
    close(nlce_file)
end

function write_to_file(
    nlce_output::AbstractDict{AbstractNLCECluster,Vector{<:Real}},
    filename::AbstractString,
)

    nlce_file = open(filename, "w")

    for (cluster, mults) in nlce_output
        write(nlce_file, "$(nv(cluster)):")
        for edge in edge_list(cluster)
            write(nlce_file, " $(join(edge, ' '))")
        end
        write(nlce_file, " : $(join(mults, ' '))\n")
    end
    close(nlce_file)
end

function write_to_file_fortran(
    nlce_output::AbstractDict{AbstractNLCECluster,Vector{<:Real}},
    filename::AbstractString,
    max_order::Integer,
)

    nlce_files = [open(filename * "_$(i).txt", "w") for i = 1:max_order]
    sorted_clusters = sort(collect(keys(nlce_output)), by = nv)

    for cluster in sorted_clusters
        edges = edge_list(cluster)
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
