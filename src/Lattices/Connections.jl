abstract type AbstractConnections end
Base.getindex(c::AbstractConnections, evs::ExpansionVertices{Int}) = _NI("getindex")

struct StrongClusterConnections <: AbstractConnections
        connections::Vector{LatticeVertices{Int}}
end

Base.getindex(c::StrongClusterConnections, evs::ExpansionVertices{Int}) = union(LatticeVertices{Int}(), c.connections[evs])

struct WeakClusterConnections <: AbstractConnections
        connections::Vector{LatticeVertices{Int}}
        rev_connections::Vector{ExpansionVertices{Int}}
        connections_matrix::Matrix{Int}
        basis_sizes::Vector{Int}
end

function Base.getindex(c::WeakClusterConnections, evs::ExpansionVertices{Int})
        lvs = union(LatticeVertices{Int}(), c.connections[evs])
        masking_matrix = .!in.(c.connections_matrix[lvs, lvs], (evs,))
        masking_matrix[diagind(masking_matrix)] .= 0
        return lvs, masking_matrix
end

just_lvs(c::WeakClusterConnections, evs::ExpansionVertices{Int}) = union(LatticeVertices{Int}(), c.connections[evs])
