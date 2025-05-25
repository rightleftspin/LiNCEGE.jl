"""
Example Ising model code from LINCEGE that reads in clusters, solves the Ising model
and writes the output to a jld file, which is useful to read into julia later
and plot.
"""

using NLCE
using JSON3
using JLD

order = 5
clusters = JSON3.read(open("examples/outputs/ex-triangular-cluster/triangular_cluster1_nn_$(order).json", "r"))
num_sites, bond_lists, multiplicities = [], [], []

for cluster in clusters
    push!(num_sites, cluster["Number of Sites"])
    push!(bond_lists, cluster["Weighted Bonds"])
    # The conversion here is nessecary because JSON3 returns arrays
    # as static by default
    push!(multiplicities, Vector(cluster["Multiplicities"]))
end

# Turn off magnetic field
B = 0
couplings = [1, 0]
temperature = range(0, 10, length = 10000)

# returns observables in 3D array, (property, temperature, order)
# In order, properties are (energy, entropy, specific heat, magnetization)
obs = NLCE.observables(
    NLCE.ising_model,
    num_sites,
    bond_lists,
    multiplicities,
    temperature,
    B,
    couplings,
    2,
    order - 1
)

# Writing all the files to the corresponding folder, creating the folder
# if it does not exist
filepath = "examples/outputs/ex-ising"
mkpath(filepath)
filename = filepath * "/triangular_cluster1_nn_obs_$(order).jld"

save(filename, "temp", temperature, "obs", obs)
