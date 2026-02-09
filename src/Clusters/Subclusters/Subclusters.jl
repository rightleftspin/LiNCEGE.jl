struct Subclusters <: AbstractClusters
    clusters::Set{ExpansionVertices}
end

function Subclusters(vs::ExpansionVertices)
    subclusters = Set{ExpansionVertices}()
    for v in vs
        push!(subclusters, ExpansionVertices(v))
    end

    Subclusters(subclusters)
end

function Subclusters(cluster::AbstractCluster, lattice::AbstractSiteExpansionLattice)
    if length(cluster) == 1
        return Subclusters()
    end

    max_order = length(cluster) - 1
    roots = Subclusters(vertices(cluster))
    visited = Subclusters()

    function try_mark(subcluster::ExpansionVertices)
        already = subcluster in visited
        if !already
           push!(visited, subcluster)
        end

        if length(subcluster) == max_order
            return false
        end
        !already
    end

    function dfs(subcluster::ExpansionVertices)
        if !try_mark(subcluster)
            return
        end

        nbrs = neighbor_subclusters(subcluster, lattice)
        for ssc in nbrs
            dfs(ssc)
        end
    end

    for sc in roots
        dfs(sc)
    end

    visited
end

_clusters(cs::Subclusters) = cs.clusters
Base.in(c::ExpansionVertices, cs::Subclusters) = c in _clusters(cs)
Base.length(cs::Subclusters) = length(_clusters(cs))
Base.iterate(cs::Subclusters) = iterate(_clusters(cs))
Base.iterate(cs::Subclusters, state) = iterate(_clusters(cs), state)
Base.push!(cs::Subclusters, c::ExpansionVertices) = push!(_clusters(cs), c)

function neighbor_subclusters(c::ExpansionVertices, lattice::AbstractSiteExpansionLattice)
    for v in c

    end
end

function find_subclusters(cs::AbstractClusters, cid::AbstractClusterID, lattice::AbstractLattice)
    subcluster_ids = Vector{AbstractClusterID}()
    for subcluster in Subclusters(cs[cid], lattice)
        # This seems particularly inefficient, may require hashing the subcluster twice, but could be saved in the state of subcluster perhaps
        if subcluster in cs
            push!(subcluster_ids, get_id(cs, subcluster))
        else
            error("Subcluster of Cluster $(cid) is not contained in the original set of clusters")
        end
    end
    subcluster_ids
end

function add_subclusters!(cs::AbstractClusters, lattice::AbstractLattice)
    for cid in cs
        add_subclusters!(cs[cid], find_subclusters(cs, cid, lattice))
    end
    cs
end
