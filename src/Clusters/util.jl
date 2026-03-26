function add_cluster!(cs::AbstractClusterSet{C,H}, c::AbstractCluster) where {C,H}
        new_cluster = C(c.evs, c.lc, ghash(cs, c))
        if new_cluster in cs
                old_cluster = pop!(cs, new_cluster)
                push!(cs, C(c.evs, c.lc + old_cluster.lc, old_cluster.ghash))
        else
                push!(cs, new_cluster)
        end
end

function clusters_from_clusters!(new_clusters::AbstractClusterSet{C,H}, old_clusters::AbstractClusterSet) where {C<:AbstractCluster,H}
        for old_cluster in old_clusters
                add_cluster!(new_clusters, old_cluster)
        end
        new_clusters
end

function clusters_from_lattice!(clusters::AbstractClusterSet{C,H}, lattice::AbstractInfiniteLattice; spawn_depth::Int=3) where {C<:AbstractCluster,H}
        max_depth = max_order(lattice)
        roots = [C(ExpansionVertices(center), clusters, lattice) for center in centers(lattice)]
        vlock = ReentrantLock()

        function try_mark(cluster::AbstractCluster)
                already = lock(vlock) do
                        in_cs = cluster in clusters
                        if !in_cs
                                push!(clusters, cluster)
                        end
                        in_cs
                end

                if length(cluster) == max_depth
                        return false
                end
                !already
        end

        function dfs(cluster::AbstractCluster, depth)
                if !try_mark(cluster)
                        return
                end

                if depth < spawn_depth
                        for ev in neighbors(lattice, cluster.evs)
                                dfs(C(union(cluster.evs, ExpansionVertices(ev)), clusters, lattice), depth + 1)
                        end
                else
                        tasks = Task[]
                        first = true
                        for ev in neighbors(lattice, cluster.evs)
                                neighbor_cluster = C(union(cluster.evs, ExpansionVertices(ev)), clusters, lattice)

                                if first
                                        dfs(neighbor_cluster, depth + 1)
                                        first = false
                                else
                                        push!(tasks, @spawn dfs(neighbor_cluster, depth + 1))
                                end
                        end
                        for t in tasks
                                fetch(t)
                        end
                end
        end

        for c in roots
                dfs(c, 0)
        end

        clusters
end
