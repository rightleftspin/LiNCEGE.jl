module NLCE

using NautyGraphs
using Graphs

# Add the basic pipeline
include("pipeline/1-grow.jl")
include("pipeline/2-prune.jl")
include("pipeline/3-combine.jl")

# Add the relevant helper functions
include("helpers/pruning.jl")

export
    isomorphic_tagging


# Below is an example of the NLCE process for a square lattice with only nearest neighbors
# up till order 4 it uses a wrapper function for speed, check the helper functions for 
# more info

# add write to file helper
#using AlgebraicNumbers: AlgebraicNumber as AN
#using Plots; gr()
#
### Tasks before next meeting
## do square, triangle and kagome lattice ising models
#
### Setting up the lattice geometry
#basis = [[AN(0), AN(0)], [AN(1), AN(0)], [AN(1)/2, sqrt(AN(3))/2]] 
#primitive_vectors = [[AN(2), AN(0)], [AN(1), sqrt(AN(3))]]
#final_order = 7
#
#for order in 1:final_order
#    println(length(simple_nlce(basis, primitive_vectors, order)))
#end

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
#    for (cluster, mult) in nlce_weights
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
