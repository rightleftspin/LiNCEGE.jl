abstract type AbstractConnections end
connections(c::AbstractConnections, evs::ExpansionVertices{Int}) = _NI("connections")

struct StrongClusterConnections <: AbstractConnections
        connections::Vector{LatticeVertices{Int}}
end

connections(c::StrongClusterConnections, evs::ExpansionVertices{Int}) = union(LatticeVertices{Int}(), c.connections[evs])
