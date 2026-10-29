# This script is automatically generated without explanatory comments.
# The documentation for this script exists at https://rightleftspin.github.io/LiNCEGE.jl/dev/examples/square_lattice.html

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

using Plots, ColorSchemes
image_unit_cell(square_uc)

image_lattice(square_uc)

m_order = 3
lattice = SiteExpansionLattice(m_order, square_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
reduce_clusters!(iso_clusters, trans_clusters)

expansion = Expansion(iso_clusters, lattice)
summation!(expansion, m_order)

json_path = joinpath(output_dir, "square_lattice_site_expansion.json")
write_to_json(expansion, lattice, json_path)

txt_path = joinpath(output_dir, "square_lattice_site_expansion.txt")
write_to_txt_file(expansion, lattice, txt_path)

temperatures = collect(range(0.5, 5.0, length=200))

clusters = import_from_json(json_path)
solvers = [IsingSolver(c) for c in clusters]

nlce_result = perform_nlce(solvers, temperatures)
resummed = apply_resummations(nlce_result, [(:Euler, 2), (:Wynn, 1)])

T_idx = argmax(temperatures .>= 2.0)
(temperature = temperatures[T_idx],
 bare = nlce_result.specific_heat[end][T_idx],
 resummed = resummed.specific_heat[end][T_idx])

using Plots

T = nlce_result.temperatures
n_orders = length(nlce_result.energy)

energy_cumsum  = hcat([sum(nlce_result.energy[1:n]) for n in 1:n_orders]...)
entropy_cumsum = hcat([sum(nlce_result.entropy[1:n]) for n in 1:n_orders]...)
cv_cumsum      = hcat([sum(nlce_result.specific_heat[1:n]) for n in 1:n_orders]...)

plot_range = max(1, n_orders - 2):n_orders
labels = reshape(["Order $n" for n in plot_range], 1, :)

p1 = plot(T, energy_cumsum[:, plot_range]; label=false, xlabel="T/J", title="Energy")
p2 = plot(T, entropy_cumsum[:, plot_range]; label=labels, xlabel="T/J", title="Entropy", legend=:top)
p3 = plot(T, cv_cumsum[:, plot_range]; label=false, xlabel="T/J", title="Specific Heat")

plot(p1, p2, p3; layout=(1, 3), size=(900, 300))
