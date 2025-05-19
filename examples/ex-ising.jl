using JLD
using NLCE

order = 6
num_sites, bond_lists, multiplicities = load("./outputs/Square_Lattice_exp_$(order).jld", "num_sites"),
load("./outputs/Square_Lattice_exp_$(order).jld", "bond_lists"),
load("./outputs/Square_Lattice_exp_$(order).jld", "multiplicities")

B = 0
couplings =  [1, 5e-2]
temperature = range(0, 10, length=100)

obs = NLCE.ising_observables(num_sites, bond_lists, multiplicities, temperature, B, couplings, 3, 4)

save("./outputs/square_lattice_exp_obs_$(order).jld", "temp", temperature, "obs", obs)
