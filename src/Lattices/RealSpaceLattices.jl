struct RealSpaceLattice{
        C<:AbstractMatrix,
        S<:AbstractMatrix{<:Integer},
        T<:AbstractVector{ExpansionVertices},
        L<:AbstractVector{<:Integer},
}
        coordinates::C
        labels::L
        translational_labels::L
        rev_connections::T
end

function RealSpaceLattice(nlce_tiling::Tiling, max_order::Int)

end

pairwise_distance(lattice::RealSpaceLattice) = pairwise(euclidean, lattice.coords, dims=1)
coords(lattice::RealSpaceLattice) = eachrow(lattice.coords)
labels(lattice::RealSpaceLattice) = lattice.labels
translational_labels(lattice::RealSpaceLattice) = lattice.translational_labels

function matching_exp_vertex(lattice::RealSpaceLattice, i::Int, j::Int)
        # Sometimes, there are dangling bonds in the lattice, however, these are usually
        # far out enough to where you will never expand into then, so they can safely
        # be set to 0 and not be dealt with
        v = intersect(lattice.rev_connections[i], lattice.rev_connections[j])
        if isempty(v)
                return zero(eltype(v))
        end

        pop!(v)
end
