struct ExpansionBond
        site1::Vector{Int}
        site2::Vector{Int}
        direction::Vector{Int}
        bond_type::Int
end

neighbor_site(bond::ExpansionBond, coordinate::AbstractVector{Int}) = [coordinate[1:end-2] + bond.direction; bond.site2[1]; bond.site2[2]]

# coordinates are written [x1 x2 x3 ...; y1 y2 y3 ...; z1 z2 z3 ...]
struct ExpansionUnitCell
        basis::Vector{Matrix{Float64}}
        primitive_vectors::Matrix{Float64}
        site_colors::Vector{Vector{Int}}
        bonds::Vector{ExpansionBond}
        expansion_bonds::Vector{Bond}
end

function ExpansionUnitCell(basis::AbstractVector{<:AbstractVector{<:AbstractVector{Float64}}}, primitive_vectors::AbstractVector{<:AbstractVector{Float64}}, bonds::AbstractVector{ExpansionBond}, expansion_bonds::AbstractVector{Bond}, site_colors::AbstractVector{<:AbstractVector{Int}})

        @assert length(site_colors) == length(basis) "Number of site colors must match number of expansion basis elements"
        @assert all(length.(primitive_vectors) .== length(primitive_vectors)) "Primitive vectors must form a square matrix"

        ExpansionUnitCell([hcat(b...) for b in basis], hcat(primitive_vectors...), site_colors, bonds, expansion_bonds)
end

basis_size(unit_cell::ExpansionUnitCell) = size.(unit_cell.basis, 2)
dimension(unit_cell::ExpansionUnitCell) = size(unit_cell.primitive_vectors, 2)
shift_unit_cell(unit_cell::ExpansionUnitCell, shift_vector::AbstractVector{Int}) = unit_cell.primitive_vectors * shift_vector[1:end-2] + unit_cell.basis[shift_vector[end-1]][:, shift_vector[end]]
shift_unit_cell(unit_cell::ExpansionUnitCell, shift_vectors::AbstractMatrix{Int}) = stack(v -> shift_unit_cell(unit_cell, v), eachcol(shift_vectors))
