
"""
Resummation techniques and related functions. These resummation
techniques will take in a multidimensional array where the first
index is the property (energy, entropy, etc), the second index
is the independent variable (temperature, magnetic field, etc),
and the third index is over each order in the NLCE bare sum,
where each order corresponds to the partial sum up till that
order, from the weights returned by LINCEGE
"""
# TODO: Write the resummation functions more generally
function bincoeff(n, k)
    r = 1
    if (k > n)
        return 0
    end

    for d = 1:k
        r = r * n / d
        n -= 1
    end

    r
end

function break_properties(properties)
    broken_props = zero(properties)
    broken_props[:, :, 1] = properties[:, :, 1]

    for i = 2:size(properties, 3)
        broken_props[:, :, i] = properties[:, :, i] - properties[:, :, i-1]
    end

    broken_props
end

function euler_resummation(properties, start)
    broken_props = break_properties(properties[:, :, 1:(end-2)])
    partial_prop = broken_props[:, :, start:end]
    out = zero(partial_prop)
    for i = 1:size(partial_prop, 3)
        delta = zero(partial_prop[:, :, 1])
        for j = 1:i
            coeff = bincoeff(i, j)
            delta += (-1)^(j) * coeff .* abs.(partial_prop[:, :, j])
        end
        out[:, :, i] += (0.5^(i + 1)) .* delta
    end

    (sum(out, dims = 3) + properties[:, :, start-1])
end

function eps_wynn(k, n, properties)
    if k == 0
        properties[:, :, n]
    elseif k == -1
        zero(properties[:, :, n])
    else
        first = eps_wynn(k - 2, n + 1, properties)
        second = eps_wynn(k - 1, n + 1, properties) .- eps_wynn(k - 1, n, properties)
        total = first + 1 ./ second
        total
    end
end

function wynn_resummation(properties, wynn_cycles)
    final_order = size(properties, 3)

    eps_wynn((2 * wynn_cycles), final_order - (2 * wynn_cycles), properties)
end
