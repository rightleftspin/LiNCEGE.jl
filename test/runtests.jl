using Test
using Plots
using JSON3
using LINCEGE
import LINCEGE: basis_size, dimension, shift_unit_cell,
        centers, max_order, n_unique_sites, neighbors, bond_matrix

# ---------------------------------------------------------------------------
# Shared unit cell definitions reused across multiple testsets
# ---------------------------------------------------------------------------

# Square Lattice Site Expansion
square_basis = [[0.0, 0.0]]
square_pvecs = [[1.0, 0.0], [0.0, 1.0]]
square_bonds = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
square_symmetries::Vector{Matrix{Float64}} = [[1 0; 0 1], [0 1; -1 0], [-1 0; 0 -1], [0 -1; 1 0], [0 1; 1 0], [0 -1; -1 0], [1 0; 0 -1], [-1 0; 0 1]]
square_uc = UnitCell(square_basis, square_pvecs, square_bonds, [1])

# Kagome Lattice Site Expansion
kagome_basis = [[0.0, 0.0], [1.0, 0.0], [0.5, sqrt(3) / 2]]
kagome_pvecs = [[2.0, 0.0], [1.0, sqrt(3)]]
kagome_bonds = [Bond(1, 2, [0, 0], 1), Bond(2, 3, [0, 0], 1), Bond(3, 1, [0, 0], 1),
        Bond(1, 2, [-1, 0], 1), Bond(1, 3, [0, -1], 1), Bond(2, 3, [1, -1], 1)]
kagome_uc = UnitCell(kagome_basis, kagome_pvecs, kagome_bonds, [1, 1, 1])

# Shastry-Sutherland Lattice Site Expansion
ss_pvecs = [[(sqrt(3) + 1) / sqrt(2), 0.0], [0.0, (sqrt(3) + 1) / sqrt(2)]]
ss_basis = [[0.0, 0.0],
        [(sqrt(3) + 1) / (2 * sqrt(2)), (sqrt(3) - 1) / (2 * sqrt(2))],
        [sqrt(2) / 2, sqrt(6) / 2],
        [(sqrt(3) + 3) / (2 * sqrt(2)), (sqrt(3) + 1) / (2 * sqrt(2))]]
ss_bonds = [Bond(1, 2, [0, 0], 1), Bond(2, 3, [0, 0], 1), Bond(1, 4, [-1, 0], 1), Bond(3, 4, [-1, 0], 1),
        Bond(3, 4, [0, 0], 2), Bond(2, 3, [0, -1], 2), Bond(1, 2, [-1, 0], 2), Bond(1, 4, [-1, -1], 2),
        Bond(2, 4, [-1, -1], 3), Bond(1, 3, [-1, 0], 3)]
ss_uc = UnitCell(ss_basis, ss_pvecs, ss_bonds, [1, 2, 3, 4])

# One Fifth Depleted Square Lattice Site Expansion
one_fifth_pvecs = [[2.0, 1.0], [-1.0, 2.0]]
one_fifth_basis = [[0.0, 0.0], [0.0, 1.0], [0.0, 2.0], [1.0, 1.0]]
one_fifth_bonds = [Bond(1, 2, [0, 0], 1), Bond(1, 3, [0, -1], 1),
        Bond(2, 4, [0, 0], 1), Bond(3, 4, [0, 1], 1),
        Bond(2, 3, [0, 0], 2), Bond(1, 4, [-1, 0], 2)]
one_fifth_uc = UnitCell(one_fifth_basis, one_fifth_pvecs, one_fifth_bonds, [1, 1, 1, 1])

# Cubic Lattice Site Expansion
cube_pvecs = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
cube_bonds = [Bond(1, 1, [1, 0, 0], 1), Bond(1, 1, [0, 1, 0], 1), Bond(1, 1, [0, 0, 1], 1)]
cube_uc = UnitCell([[0.0, 0.0, 0.0]], cube_pvecs, cube_bonds, [1])




# ---------------------------------------------------------------------------

@testset verbose = true "LINCEGE.jl" begin
        include("test_vertices.jl")
        include("test_unitcells.jl")
        include("test_lattices.jl")
        include("test_clusters.jl")
        include("test_expansions.jl")
end
