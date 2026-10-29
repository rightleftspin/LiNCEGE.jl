# This script is automatically generated without explanatory comments.
# The documentation for this script exists at https://rightleftspin.github.io/LiNCEGE.jl/dev/examples/one_fifth_depleted_lattice.html

output_dir = joinpath(@__DIR__, "outputs")
mkpath(output_dir)

using Lincege
if !@isdefined(output_dir)
    output_dir = joinpath(pkgdir(Lincege), "examples", "outputs")
end
mkpath(output_dir)

using Lincege

one_fifth_basis = [[0.0, 0.0], [0.0, 1.0], [0.0, 2.0],[1.0,1.0]]
one_fifth_pvecs = [[2.0, 1.0], [-1.0, 2.0]]
one_fifth_bonds = [Bond(1,2,[0,0],1),Bond(1,3,[0,-1],1),Bond(2,4,[0,0],1), Bond(3,4,[0,1],1),Bond(2,3,[0,0],2),Bond(1,4,[-1,0],2)]
one_fifth_uc = UnitCell(one_fifth_basis, one_fifth_pvecs, one_fifth_bonds, [1, 1, 1, 1])

using Plots, ColorSchemes
image_unit_cell(one_fifth_uc)

image_lattice(one_fifth_uc)

m_order = 3
lattice = SiteExpansionLattice(m_order, one_fifth_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
reduce_clusters!(iso_clusters, trans_clusters)

expansion = Expansion(iso_clusters, lattice)
summation!(expansion, m_order)

write_to_json(expansion, lattice, joinpath(output_dir, "one_fifth_depleted_lattice_site_expansion.json"))
