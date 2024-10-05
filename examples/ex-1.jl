module ex1
using NLCE
"""
Example generating a square lattice and writing to a fortran formatted file
"""
basis = [[0, 0]]
primitive_vec = [[1, 0], [0, 1]]
neighborhood = [1]
max_order = 4

nlce_clusters = simple_NLCE(basis, primitive_vec, neighborhood, max_order)

filename = "examples/outputs/square_nn"
write_to_file_fortran(nlce_clusters, filename, max_order)

#xmax = 2
## Initialize the temperature grid
#total_temp = collect(0.05:0.05:xmax)
#all_energies = []
#
## Loop over each desired order
#for order in 5:final_order
#    # Initialize the energies
#    total_energy = zeros(length(total_temp))
#    # Perform the NLCE sum and find resulting weights
#    nlce_weights = simple_nlce(basis, primitive_vectors, order)
#    for (cluster, mult) in ( t = findall(p, x); 
#        # Find the energies for each cluster and add it to the total sum
#        ising_energy = ising_energies(cluster)
#        total_energy += mult .* energy_solver(total_temp, ising_energy) 
#    end
#    push!(all_energies, total_energy)
#end
#
## Plot
#plot(total_temp, all_energies, label = collect(7:final_order))
#xlims!(0, xmax)
#ylims!(-1, 0)
#title!("Ising Energy for a triangular Lattice")
#xlabel!("Temperature")
#ylabel!("Energy")
#savefig("ising_energy.png")
#
end
