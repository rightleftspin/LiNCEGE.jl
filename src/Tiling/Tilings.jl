struct Tiling{
        B<:AbstractVector{<:AbstractVector{<:AbstractVector{<:Real}}},
        T<:AbstractVector{AbstractVector{<:Real}},
        L<:AbstractVector{AbstractVector{<:Integer}},
        N<:AbstractVector{<:Real},
        F,
}
        tiling_unit::B
        translation_vectors::T
        translation_labels::L
        expansion_neighbors::N
        real_space_neighbors::N
        labels::L
        neighbor_fn::F
end

function Tiling(
        basis::AbstractVector{<:AbstractVector{<:Real}},
        primitive_vectors::AbstractVector{AbstractVector{<:Real}},
        neighbors::AbstractVector{<:Real};
        labels::AbstractVector{<:Integer}=repeat([1], length(basis)),
        neighbor_fn::Function=(
                (
                        lattice::RealSpaceLattice,
                        adj_coordinates::Vector{CartesianIndex{2}},
                        distance_index::Int,
                ) -> repeat([distance_index], length(adj_coordinates))
        ),
        expand_by_basis::Bool=false,
)

        tiling_unit = Vector{Vector{Vector{Float64}}}[]
        if expand_by_basis
                tiling_unit = [basis]
        else
                for basis_elem in basis
                        push!(tiling_unit, [basis_elem])
                end
        end

        translation_vectors = primitive_vectors
        translation_labels = [collect(1:length(basis))]
        opt_labels = [labels]

        Tiling(
                tiling_unit,
                translation_vectors,
                translation_labels,
                neighbors,
                neighbors,
                opt_labels,
                neighbor_fn,
        )
end

dimension(tiling::Tiling) = length(tiling.translation_vectors[1])
real_space_neighbors(tiling::Tiling) = tiling.real_space_neighbors
