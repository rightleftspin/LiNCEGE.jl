function get_subgraphs(c::AbstractCluster, lattice::AbstractLattice)
        if length(c) == 1
                return Set()
        end

        max_depth = length(c) - 1
        ctrs = c.vs
        roots = [typeof(ctrs)(center) for center in ctrs]
        visited = Set()

        function try_mark(cluster)
                already = cluster in visited
                if !already
                        push!(visited, cluster)
                end

                if length(cluster) == max_depth
                        return false
                end
                !already
        end

        function dfs(cluster)
                if !try_mark(cluster)
                        return
                end
                for v in neighbors(lattice, cluster)
                        if v in c.vs
                                dfs(union(cluster, typeof(ctrs)(v)))
                        end
                end
        end

        for c in roots
                dfs(c)
        end

        visited
end

function adj_mat_to_edge_list(adj_matrix::AbstractMatrix{<:Real})
        edge_list = []

        n = size(adj_matrix, 1)
        for i in 1:n
                for j in i+1:n
                        weight = adj_matrix[i, j]
                        if weight != 0
                                push!(edge_list, (i, j, weight))
                        end
                end
        end

        edge_list
end

"""
    write_to_json(expansion, lattice, filepath)

Serialize an `Expansion` and its associated `lattice` geometry to a JSON file at
`filepath`.  The file contains a JSON array — one object per cluster — with the
following fields:

- `cluster_hash`   — unique cluster identifier (string representation of UInt64)
- `order`          — 1-based order index into the expansion's `order_ids`
- `n_sites`        — number of sites in the cluster
- `lattice_constant` — lattice constant used in the NLCE summation
- `coordinates`    — list of Cartesian coordinate vectors (one per site)
- `site_colors`    — list of integer site-color labels (one per site)
- `bonds`          — edge list `[i, j, weight]` using local 1-based indices
- `subgraphs`      — list of subgraph hashes (string representation of UInt64)
- `weights`        — dict mapping cluster hash strings to NLCE weight floats
"""
function write_to_json(e::Expansion, lattice::AbstractLattice, filepath::String)
        all_coords = get_coordinates(lattice)
        all_colors = get_site_colors(lattice)
        adj = bond_matrix(lattice)

        clusters_data = []

        for (order_idx, cluster_hashes) in enumerate(e.order_ids)
                for cluster_hash in cluster_hashes
                        cluster = e.expansion_clusters[cluster_hash]
                        vs = cluster.vertices
                        n = length(vs)

                        coords = [collect(col) for col in eachcol(all_coords[:, vs])]
                        colors = collect(all_colors[vs])
                        bonds = [[b[1], b[2], b[3]] for b in adj_mat_to_edge_list(adj[vs, vs])]
                        sgs = [string(sg) for sg in cluster.subgraphs]
                        wts = Dict(string(k) => v for (k, v) in cluster.weights)

                        push!(clusters_data, Dict(
                                "cluster_hash" => string(cluster_hash),
                                "order" => order_idx,
                                "n_sites" => n,
                                "lattice_constant" => cluster.lattice_constant,
                                "coordinates" => coords,
                                "site_colors" => colors,
                                "bonds" => bonds,
                                "subgraphs" => sgs,
                                "weights" => wts
                        ))
                end
        end

        open(filepath, "w") do io
                JSON.print(io, clusters_data, 2)
        end
end