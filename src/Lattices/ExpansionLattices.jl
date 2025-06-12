struct ExpansionLattice{V<:AbstractVector{ExpansionVertices},C<:AbstractVector{RealSpaceVertices}}
        adjacency_list::V
        connections::C
end

function ExpansionLattice(nlce_tiling::Tiling, max_order::Int)

end

neighbors(lattice::ExpansionLattice, vertex::Integer) = lattice.adjacency_list[vertex]

real_space_vertices(lattice::ExpansionLattice, vertex::Unsigned) =
        lattice.connections[vertex]
real_space_vertices(lattice::ExpansionLattice, exp_v::ExpansionVertices) =
        sorted_real_space_vertices(lattice::ExpansionLattice, vertex::Unsigned) =
                lattice.connections[vertex]
sorted_real_space_vertices(lattice::ExpansionLattice, exp_v::ExpansionVertices) =
        union(lattice.connections[exp_v])

get_mask(lattice::ExpansionLattice, g::ExpansionLatticeGraph, vertex::Unsigned) =
        get_mask(g, vertex, real_space_vertices(lattice, vertex))
get_mask(
        lattice::ExpansionLattice,
        g::ExpansionLatticeGraph,
        exp_v::ExpansionVertices,
) = get_mask(g, exp_v, real_space_vertices(lattice, exp_v))
