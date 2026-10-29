# This script is automatically generated without explanatory comments.
# The documentation for this script exists at https://rightleftspin.github.io/LiNCEGE.jl/dev/examples/finite_square_lattice.html

output_dir = joinpath(@__DIR__, "outputs")
mkpath(output_dir)

using Lincege
if !@isdefined(output_dir)
    output_dir = joinpath(pkgdir(Lincege), "examples", "outputs")
end
mkpath(output_dir)

using Lincege

square_basis = [[0.0, 0.0]]
square_pvecs = [[1.0, 0.0], [0.0, 1.0]]
square_bonds = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
square_uc = UnitCell(square_basis, square_pvecs, square_bonds, [1])

finite_order = 3
finite_lattice = FiniteLattice((4, 4), square_uc, finite_order; periodic=true)

finite_connected = ConnectedClusterSet(finite_lattice)
clusters_from_lattice!(finite_connected, finite_lattice)

finite_symmetric = SymmetricClusterSet(finite_lattice, :Square)
reduce_clusters!(finite_symmetric, finite_connected)

finite_isomorphic = IsomorphicClusterSet(finite_lattice)
reduce_clusters!(finite_isomorphic, finite_connected)

finite_expansion = Expansion(finite_isomorphic, finite_lattice)
summation!(finite_expansion, finite_order)

finite_json_path = joinpath(output_dir, "square_lattice_periodic_finite_cluster_expansion.json")
write_to_json(finite_expansion, finite_lattice, finite_json_path)
