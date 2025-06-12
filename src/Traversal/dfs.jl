struct Consensus
        is_done::AbstractVector
        thread_id::Int
end

all_threads_done(c::Consensus) = (sum(c.is_done) == 0)
thread_empty!(c::Consensus) = (c.is_done[c.thread_id] = 0)
thread_working!(c::Consensus) = (c.is_done[c.thread_id] = 1)


function dfs_kernel!(
        dfs_lock::ReentrantLock,
        next::Stack{AbstractCluster},
        lattice::ExpansionLattice,
        parents::AbstractSet{Subgraph},
        max_order::Integer,
        c::Consensus,
)

        while true

                lock(dfs_lock)
                if all_threads_done(c)
                        unlock(dfs_lock)
                        break
                end

                if isempty(next)
                        println("thread $(c.thread_id) empty")
                        thread_empty!(c)
                        println("$(c.is_done) state")
                        unlock(dfs_lock)
                        sleep(0.1)
                        continue
                end

                thread_working!(c)
                subgraph = pop!(next)
                println("thread $(c.thread_id) working")

                unlock(dfs_lock)
                if length(subgraph) < max_order
                        filtered_subgraphs = filter(!in(parents), neighbor_subgraphs(subgraph, lattice))
                        lock(dfs_lock)
                        for fs in filtered_subgraphs
                                push!(parents, fs)
                                push!(next, fs)
                        end
                        unlock(dfs_lock)

                end
        end
        println("Finished")
end

function dfs!(
        next::Stack{Subgraph},
        cluster::Cluster,
        source::Subgraph,
        parents::AbstractSet{Subgraph},
        max_order::Integer,
)

        push!(next, source)
        push!(parents, source)

        dfs_lock = ReentrantLock()

        con_array = ones(Int, Threads.nthreads())

        Threads.@threads for i = 1:Threads.nthreads()
                c = Consensus(view(con_array, :), Threads.threadid())
                dfs_kernel!(dfs_lock, next, cluster, parents, max_order, c)
        end

        parents
end

function grow_lower_site!(
        cluster::Cluster,
        parents::AbstractSet{Subgraph},
        start::Integer,
        max_order::Integer,
)

        root_subgraph = Subgraph(start, cluster)

        next = Stack{Subgraph}()

        dfs!(next, cluster, root_subgraph, parents, max_order)
end

function grow_lower(cluster::Cluster, start::Vector{<:Integer}, max_order::Integer)

        parents = Set{Subgraph}()

        for vertex in start
                grow_lower_site!(cluster, parents, vertex, max_order)
        end

        parents
end
