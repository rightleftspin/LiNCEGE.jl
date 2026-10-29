# This script is automatically generated without explanatory comments.
# The documentation for this script exists at https://rightleftspin.github.io/LiNCEGE.jl/dev/examples/shastry-sutherland.html

output_dir = joinpath(@__DIR__, "outputs")
mkpath(output_dir)

using Lincege
if !@isdefined(output_dir)
    output_dir = joinpath(pkgdir(Lincege), "examples", "outputs")
end
mkpath(output_dir)

using Lincege

shsu_basis = [[0.0, 0.0], [(sqrt(3) + 1) / (2 * sqrt(2)), (sqrt(3) - 1) / (2 * sqrt(2))], [sqrt(2) / 2, sqrt(6) / 2], [(sqrt(3) + 3) / (2 * sqrt(2)), (sqrt(3) + 1) / (2 * sqrt(2))]]
shsu_pvecs = [[(sqrt(3) + 1) / sqrt(2), 0.0], [0.0, (sqrt(3) + 1) / sqrt(2)]]
shsu_bonds = [Bond(1, 2, [0, 0], 1), Bond(2, 3, [0, 0], 1), Bond(1, 4, [-1, 0], 1), Bond(3, 4, [-1, 0], 1),
        Bond(3, 4, [0, 0], 2), Bond(2, 3, [0, -1], 2), Bond(1, 2, [-1, 0], 2), Bond(1, 4, [-1, -1], 2),
        Bond(2, 4, [-1, -1], 3), Bond(1, 3, [-1, 0], 3)]
shsu_uc = UnitCell(shsu_basis, shsu_pvecs, shsu_bonds, [1, 2, 3, 4])

using Plots, ColorSchemes
image_unit_cell(shsu_uc)

image_lattice(shsu_uc)

m_order = 3
lattice = SiteExpansionLattice(m_order, shsu_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
reduce_clusters!(iso_clusters, trans_clusters)

expansion = Expansion(iso_clusters, lattice)
summation!(expansion, m_order)

write_to_json(expansion, lattice, joinpath(output_dir, "shastry_sutherland_lattice_site_expansion.json"))
