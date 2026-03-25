function generate_cartesian_coordinates(dimension::Int, half_side_length::Int)
        # Forces the lattice to have a strict center point
        diameter = 2 * half_side_length + 1
        # Total number of coordinates for the entire lattice
        max_coords = diameter^dimension
        coords = repeat(0:(max_coords-1), 1, dimension)'

        for dim = 0:(dimension-1)
                coords[dim+1, :] =
                        div.(coords[dim+1, :], diameter^dim) .% diameter .- half_side_length
        end

        coords
end

function generate_coordinates(max_order::Int, num_basis_elements::Int, dimension::Int)
        primitive_coordinates = generate_cartesian_coordinates(dimension, max_order)

        hcat([vcat(primitive_coordinates, repeat([i], size(primitive_coordinates, 2))') for i in 1:num_basis_elements]...)
end

function generate_lattice_data(coordinates::Matrix{Int}, unit_cell::UnitCell)
        n = size(coordinates, 2)

        coord_index = Dict{Vector{Int},Int}()
        for (i, col) in enumerate(eachcol(coordinates))
                coord_index[collect(col)] = i
        end

        adj_matrix = zeros(Int, n, n)
        neighbor_sets = [BitSet() for _ in 1:n]

        for (index, col) in enumerate(eachcol(coordinates))
                coord_vec = collect(col)
                for (neighbor_coord, bond) in zip(find_possible_neighbors(unit_cell, coord_vec), unit_cell.bonds)
                        ni = get(coord_index, neighbor_coord, 0)
                        if ni != 0
                                adj_matrix[index, ni] = bond.bond_type
                                adj_matrix[ni, index] = bond.bond_type
                                push!(neighbor_sets[index], ni)
                                push!(neighbor_sets[ni], index)
                        end
                end
        end

        neighbor_list = [ExpansionVertices{Int}(bs) for bs in neighbor_sets]
        return adj_matrix, neighbor_list
end

function find_centers(coordinates::AbstractMatrix{Int})
        center_indices = findall(col -> all(x -> x == 0, col[1:end-1]), eachcol(coordinates))
        ExpansionVertices(center_indices)
end
