
function ising_energies(cluster)
    cluster_size = nv(cluster)
    energies = Vector{Float64}(undef, (2^cluster_size))
    for i = 0:((2^cluster_size)-1)
        sites = digits(i, base = 2, pad = cluster_size)
        energy = 0
        for bond in edges(cluster)
            if src(bond) < dst(bond)
                if (sites[src(bond)] == sites[dst(bond)])
                    energy += 0.25
                else
                    energy -= 0.25
                end
            end
        end
        energies[i+1] = energy
    end
    return energies
end

function energy_solver(temps, energies)
    temps_en = copy(collect(temps))
    #    temps_free_en = copy(collect(temps))
    #    temps_en_sq = copy(collect(temps))
    for (index, temp) in enumerate(temps)
        partition = exp.(-energies ./ temp)
        temps_en[index] = sum(energies .* partition) / sum(partition)
        #        temps_free_en[index] = sum((log.(partition)) .* partition) / sum(partition)
        #        temps_en_sq[index] = sum((energies.^2) .* partition) / sum(partition)
    end

    #return temps_en, temps_free_en, temps_en_sq
    return temps_en
end
