
struct RealSpaceLattice{C<:AbstractMatrix,S<:AbstractMatrix{<:Integer}}
        coordinates::C
        sublattice_coordinates::S
end

struct ExpansionLattice{T<:AbstractVector{<:AbstractVector{<:Integer}},M<:AbstractMatrix{<:Integer}}
        adjacency_list::T
        connections::T
        expansion_lattice_graph::ExpansionLatticeGraph{M}
end

neighbors(lattice::ExpansionLattice, vertex::Integer) = lattice.adjacency_list[vertex]
real_space_vertices(lattice::ExpansionLattice, vertex::Integer) = lattice.connections[vertex]
real_space_vertices(lattice::ExpansionLattice, vertices::AbstractSet{<:Integer}) = lattice.connections[vertices]
get_mask(lattice::ExpansionLattice, vertex::Integer) = get_mask(lattice.expansion_lattice_graph, vertex, real_space_vertices(vertex))
get_mask(lattice::ExpansionLattice, vertices::AbstractSet{<:Integer}) = get_mask(lattice.expansion_lattice_graph, vertices, real_space_vertices(vertices))
