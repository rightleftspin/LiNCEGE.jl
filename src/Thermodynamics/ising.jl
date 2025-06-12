"""
Utility functions for testing various models
"""
function ising_model(num_sites, bonds, B, couplings)
    J, mu = couplings
    num_spins = 2 ^ num_sites
    energies = zeros(num_spins)
    magnetizations = zeros(num_spins)

    for spin_config = 0:(num_spins-1)
        spins = (2 * digits(spin_config, base = 2, pad = num_sites)) .- 1
        energies[spin_config+1] -= B * mu * sum(spins)
        for bond in bonds
            if spins[bond[1]] == spins[bond[2]]
                energies[spin_config+1] += J
            else
                energies[spin_config+1] -= J
            end
        end
        magnetizations[spin_config+1] = sum(spins)
    end

    energies, magnetizations
end

function get_bit(value, index)
    (value >> index) & 1
end

function toggle_bit(value, index)
    value ⊻ (1 << index)
end

function toggle_bits(value, index1, index2)
    (value ⊻ (1 << index1)) ⊻ (1 << index2)
end

function get_bits(value, num_sites)
    ([get_bit(value, i) for i = 0:(num_sites-1)])
end

# For a spin 1/2 system, this requires the couplings, along with the crystal field splitting
# the number of sites in the cluster, and bonds, returns the hamiltonian matrix
function xyz_hamil(couplings, d, g, B, mu_b, number_sites, bonds)
    j11, j22, j33 = couplings[1, 1], couplings[2, 2], couplings[3, 3]
    j12, j23, j13 = couplings[1, 2], couplings[2, 3], couplings[1, 3]
    j21, j32, j31 = couplings[2, 1], couplings[3, 2], couplings[3, 1]

    states = 0:(2^number_sites-1)

    hamil = zeros(Complex, 2 ^ number_sites, 2 ^ number_sites)

    for state in states
        # Site Interactions
        for (ind, site) in enumerate(get_bits(state, number_sites))
            spin = (2 * site) - 1
            # Magnetic Field
            hamil[state+1, toggle_bit(state, ind-1)+1] -= g * mu_b * B * 0.5
            # Crystal Field Splitting
            hamil[state+1, state+1] += spin * d * 0.5
        end

        # Bond Interactions
        for bond in bonds
            site_1 = get_bit(state, bond[1] - 1)
            site_2 = get_bit(state, bond[2] - 1)

            if site_1 == 1
                if site_2 == 1
                    hamil[state+1, toggle_bits(state, bond[1]-1, bond[2]-1)+1] += j11 / 4
                    hamil[state+1, toggle_bits(state, bond[1]-1, bond[2]-1)+1] -= j22 / 4
                    hamil[state+1, state+1] += j33 / 4

                    hamil[state+1, toggle_bits(state, bond[1]-1, bond[2]-1)+1] +=
                        j12 / (4im)
                    hamil[state+1, toggle_bits(state, bond[1]-1, bond[2]-1)+1] +=
                        j21 / (4im)
                    hamil[state+1, toggle_bit(state, bond[1]-1)+1] -= j13 / 4
                    hamil[state+1, toggle_bit(state, bond[2]-1)+1] -= j31 / 4
                    hamil[state+1, toggle_bit(state, bond[1]-1)+1] -= j23 / (4im)
                    hamil[state+1, toggle_bit(state, bond[2]-1)+1] -= j32 / (4im)
                end

                if site_2 == 0
                    hamil[state+1, toggle_bits(state, bond[1]-1, bond[2]-1)+1] += j11 / 4
                    hamil[state+1, toggle_bits(state, bond[1]-1, bond[2]-1)+1] += j22 / 4
                    hamil[state+1, state+1] -= j33 / 4

                    hamil[state+1, toggle_bits(state, bond[1]-1, bond[2]-1)+1] -=
                        j12 / (4im)
                    hamil[state+1, toggle_bits(state, bond[1]-1, bond[2]-1)+1] +=
                        j21 / (4im)
                    hamil[state+1, toggle_bit(state, bond[1]-1)+1] += j13 / 4
                    hamil[state+1, toggle_bit(state, bond[2]-1)+1] -= j31 / 4
                    hamil[state+1, toggle_bit(state, bond[1]-1)+1] += j23 / (4im)
                    hamil[state+1, toggle_bit(state, bond[2]-1)+1] += j32 / (4im)
                end

            elseif site_1 == 0
                if site_2 == 1
                    hamil[state+1, toggle_bits(state, bond[1]-1, bond[2]-1)+1] += j11 / 4
                    hamil[state+1, toggle_bits(state, bond[1]-1, bond[2]-1)+1] += j22 / 4
                    hamil[state+1, state+1] -= j33 / 4

                    hamil[state+1, toggle_bits(state, bond[1]-1, bond[2]-1)+1] +=
                        j12 / (4im)
                    hamil[state+1, toggle_bits(state, bond[1]-1, bond[2]-1)+1] -=
                        j21 / (4im)
                    hamil[state+1, toggle_bit(state, bond[1]-1)+1] -= j13 / 4
                    hamil[state+1, toggle_bit(state, bond[2]-1)+1] += j31 / 4
                    hamil[state+1, toggle_bit(state, bond[1]-1)+1] += j23 / (4im)
                    hamil[state+1, toggle_bit(state, bond[2]-1)+1] += j32 / (4im)

                end
                if site_2 == 0
                    hamil[state+1, toggle_bits(state, bond[1]-1, bond[2]-1)+1] += j11 / 4
                    hamil[state+1, toggle_bits(state, bond[1]-1, bond[2]-1)+1] -= j22 / 4
                    hamil[state+1, state+1] += j33 / 4

                    hamil[state+1, toggle_bits(state, bond[1]-1, bond[2]-1)+1] -=
                        j12 / (4im)
                    hamil[state+1, toggle_bits(state, bond[1]-1, bond[2]-1)+1] -=
                        j21 / (4im)
                    hamil[state+1, toggle_bit(state, bond[1]-1)+1] += j13 / 4
                    hamil[state+1, toggle_bit(state, bond[2]-1)+1] += j31 / 4
                    hamil[state+1, toggle_bit(state, bond[1]-1)+1] -= j23 / (4im)
                    hamil[state+1, toggle_bit(state, bond[2]-1)+1] -= j32 / (4im)
                end
            end
        end
    end
    (Hermitian(hamil))
end

function xyz_eigs(num_sites, bonds, B, couplings)

    couplings_mat, d, g, mu_b = couplings

    xyz = xyz_hamil(couplings_mat, d, g, B, mu_b, num_sites, bonds)

    energies = eigvals(xyz)
    # TODO: get magnetization working
    magnetizations = zero(energies)

    (energies, magnetizations)
end

function observables(
    model_function,
    sites_per_cluster,
    bond_lists,
    multiplicities,
    temperatures,
    B,
    couplings,
    wynn_cycles,
    euler_start,
)
    resums = 2
    # Energy, Entropy, Specific Heat, Magnetization
    # Extra resum multiplicities for euler and wynn
    # Properties: property, order, temperature
    properties = zeros(4, length(temperatures), length(multiplicities[1]) + resums)

    i = 1
    for (num_sites, bond_list, mult) in zip(sites_per_cluster, bond_lists, multiplicities)
        println(i // length(sites_per_cluster))
        i += 1
        eigs, mags = model_function(num_sites, bond_list, B, couplings)

        exp_energy_temp_matrix = exp.(-transpose(eigs) ./ temperatures)
        partition_function = sum(exp_energy_temp_matrix, dims = 2)
        avg_energy = (exp_energy_temp_matrix * eigs) ./ partition_function

        mult_plus_resum = transpose(append!(mult, repeat([0], resums)))
        properties[1, :, :] += mult_plus_resum .* avg_energy
        properties[2, :, :] +=
            mult_plus_resum .* (log.(partition_function) + (avg_energy ./ temperatures))
        properties[3, :, :] +=
            mult_plus_resum .* (
                (
                    ((exp_energy_temp_matrix * (eigs .^ 2)) ./ partition_function) -
                    (avg_energy .^ 2)
                ) ./ (temperatures .^ 2)
            )
        properties[4, :, :] +=
            mult_plus_resum .* ((exp_energy_temp_matrix * mags) ./ partition_function)

    end

    properties[:, :, end-1] = euler_resummation(properties, euler_start)
    properties[:, :, end] = wynn_resummation(properties, wynn_cycles)
    properties
end
