# This script is automatically generated without explanatory comments.
# The documentation for this script exists at https://rightleftspin.github.io/LiNCEGE.jl/dev/examples/lieb_lattice.html

output_dir = joinpath(@__DIR__, "outputs")
mkpath(output_dir)

using Lincege
if !@isdefined(output_dir)
    output_dir = joinpath(pkgdir(Lincege), "examples", "outputs")
end
mkpath(output_dir)

using Lincege

lieb_basis = [[0.0, 0.0], [0.0, 1.0], [1, 0]]
lieb_pvecs = [[2.0, 0.0], [0.0, 2.0]]
lieb_bonds = [Bond(1,3,[1,0],1), Bond(1,2,[0,1],2)]
lieb_uc = UnitCell(lieb_basis, lieb_pvecs, lieb_bonds, [1, 1, 1])

using Plots, ColorSchemes
image_unit_cell(lieb_uc)

image_lattice(lieb_uc)

m_order = 3
lattice = SiteExpansionLattice(m_order, lieb_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
reduce_clusters!(iso_clusters, trans_clusters)

expansion = Expansion(iso_clusters, lattice)
summation!(expansion, m_order)

write_to_json(expansion, lattice, joinpath(output_dir, "lieb_lattice_site_expansion.json"))
