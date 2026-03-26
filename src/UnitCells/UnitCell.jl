struct Bond
        site1::Int
        site2::Int
        direction::Vector{Int}
        bond_type::Int
end

neighbor_site(bond::Bond, coordinate::AbstractVector{Int}) = [coordinate[1:end-1] + bond.direction; bond.site2]

# coordinates are written [x1 x2 x3 ...; y1 y2 y3 ...; z1 z2 z3 ...]
struct UnitCell
        basis::Matrix{Float64}
        primitive_vectors::Matrix{Float64}
        site_colors::Vector{Int}
        bonds::Vector{Bond}
end

function UnitCell(basis::AbstractVector{<:AbstractVector{Float64}}, primitive_vectors::AbstractVector{<:AbstractVector{Float64}}, bonds::AbstractVector{Bond}, site_colors::AbstractVector{Int})

        @assert length(site_colors) == length(basis) "Number of site colors must match number of sites in the basis."
        @assert all(length.(primitive_vectors) .== length(primitive_vectors)) "Primitive vectors must form a square matrix"
        @assert all(length.(basis) .== length(basis[1])) "Basis vectors must all have the same dimension"


        UnitCell(hcat(basis...), hcat(primitive_vectors...), site_colors, bonds)
end

basis_size(unit_cell::UnitCell) = size(unit_cell.basis, 2)
dimension(unit_cell::UnitCell) = size(unit_cell.primitive_vectors, 2)
shift_unit_cell(unit_cell::UnitCell, shift_vector::AbstractVector{Int}) = unit_cell.primitive_vectors * shift_vector[1:end-1] + unit_cell.basis[:, shift_vector[end]]
shift_unit_cell(unit_cell::UnitCell, shift_vectors::AbstractMatrix{Int}) = stack(v -> shift_unit_cell(unit_cell, v), eachcol(shift_vectors))
find_possible_neighbors(unit_cell::UnitCell, coordinate::AbstractVector{Int}) = [neighbor_site(bond, coordinate) for bond in unit_cell.bonds if bond.site1 == coordinate[end]]
