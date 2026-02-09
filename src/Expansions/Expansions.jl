module Expansions

import DataStructures.Deque

import LINCEGE:
    Vertices.AbstractVertices,
    Lattices.AbstractLattice,
    Clusters.AbstractClusters

abstract type AbstractExpansion end

Base.getindex(e::AbstractExpansion, i::Int, order::Int) = _NI("Base.getindex")
Base.setindex!(e::AbstractExpansion, l::Rational, i::Int, order::Int) = _NI("Base.setindex!")

function update_lattice_constant!(e::AbstractExpansion, i::Int, order::Int, lattice_constant::Rational)
    e[i, order] += lattice_constant
end

function _summation!(e::AbstractExpansion, cs::AbstractClusters, max_order::Int)
    for order in 1:max_order
        for cluster_id in get_order(cs, order)
            updated_lattice_constant = lattice_constant(cs, cluster_id)
            update_lattice_constant!(e, cluster_id, order, updated_lattice_constant)

            subgraphs = Deque{}()
            push!(subgraphs, (cluster_id, 0))
            while !(length(subgraphs) == 0)
                current_cluster_id, current_depth = pop!(subgraphs)
                for subcluster_id in get_subclusters(cs, current_cluster_id)
                    update_lattice_constant!(e, subcluster_id, order, ((-1)^(current_depth + 1)) * updated_lattice_constant)
                    subgraphs.append((subcluster_id, current_depth + 1))
                end
            end
        end
    end
    e
end

end
