"""
Example Heisenberg model code from LINCEGE that reads in clusters, solves the Ising model
and writes the output to a jld file, which is useful to read into julia later
and plot.
"""

using NLCE
using JSON3
using JLD

order = 5
clusters = JSON3.read(open("examples/outputs/ex-triangular-cluster/triangular_cluster1_nn_$(order).json", "r"))
#clusters = JSON3.read(open("examples/outputs/ex-triangular/triangular_nn_$(order).json", "r"))
num_sites, bond_lists, multiplicities = [], [], []

# Writing all the files to the corresponding folder, creating the folder
# if it does not exist
filepath = "examples/outputs/ex-heisenberg"
mkpath(filepath)
filename = filepath * "/triangular_cluster1_nn_heisenberg_obs_$(order).jld"

for cluster in clusters
    push!(num_sites, cluster["Number of Sites"])
    push!(bond_lists, cluster["Weighted Bonds"])
    # The conversion here is nessecary because JSON3 returns arrays
    # as static by default
    push!(multiplicities, Vector(cluster["Multiplicities"]))
end

# Turn off magnetic field
B = 0
# all couplings in matrix, d, g, mu_b
# since J11 = J22 = J33, this is the heisenberg model
couplings = [[1 0 0; 0 1 0; 0 0 1;], 0, 0, 0]
temperature = range(0, 5, length = 1000)

# returns observables in 3D array, (property, temperature, order)
# In order, properties are (energy, entropy, specific heat, magnetization)
obs = NLCE.observables(
    NLCE.xyz_eigs,
    num_sites,
    bond_lists,
    multiplicities,
    temperature,
    B,
    couplings,
    fld(order, 2),
    order - 1,
)

save(filename, "temp", temperature, "obs", obs)
