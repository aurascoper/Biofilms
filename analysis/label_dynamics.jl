module LabelDynamics

"""Sampled identity history on a fixed domain. No numeric label averaging.

Returns are arrivals at the MCS-0 reference following one or more recorded
different labels. End-of-window absences are incomplete episodes, not returns.
The API also accepts synthetic lattices so controls use the production path.
"""
function sampled_history(labels, mask, initial_species; stride::Int = 1)
    length(labels) >= 2 || error("at least two observations are required")
    stride >= 1 || error("stride must be positive")
    (length(labels) - 1) % stride == 0 || error("cadence must include the final observation")
    all(size(a) == size(mask) for a in labels) || error("grid mismatch")
    size(initial_species) == size(mask) || error("stratum grid mismatch")
    any(mask) || error("interior mask is empty")
    all(all(>=(0), a[mask]) for a in labels) || error("negative interior label")
    indices = collect(1:stride:length(labels))
    reference = labels[1]
    transitions = zeros(Int32, size(mask))
    returns = zeros(Int32, size(mask))
    hidden = zeros(Int32, size(mask))
    persistent = copy(mask)
    groups = [("all_interior", -1, copy(mask))]
    for s in 0:7
        push!(groups, ("initial_species", s, mask .& (initial_species .== s)))
    end
    for label in sort(unique(reference[mask]))
        push!(groups, ("initial_label", Int(label), mask .& (reference .== label)))
    end
    curves = [Int[count(g[3])] for g in groups]
    for (j, current) in enumerate(indices[2:end])
        previous = indices[j]
        changed = mask .& (labels[current] .!= labels[previous])
        transitions .+= changed
        # previous != reference also proves a sampled departure exists.
        returns .+= mask .& (labels[previous] .!= reference) .& (labels[current] .== reference)
        persistent .&= labels[current] .== reference
        any_fine_change = falses(size(mask))
        for fine in (previous + 1):current
            any_fine_change .|= mask .& (labels[fine] .!= labels[fine - 1])
        end
        hidden .+= any_fine_change .& (labels[current] .== labels[previous])
        for (k, (_, _, support)) in enumerate(groups)
            push!(curves[k], count(persistent .& support))
        end
    end
    rows = [Dict("stratum" => kind, "initial_label" => label,
                 "denominator_sites" => count(support),
                 "transition_count" => sum(transitions[support]),
                 "return_episode_count" => sum(returns[support]),
                 "sites_with_transition" => count(>(0), transitions[support]),
                 "sites_with_return" => count(>(0), returns[support]),
                 "hidden_reversal_intervals" => sum(hidden[support]),
                 "persistent_site_counts" => curves[k],
                 "persistence_fraction" => [count(support) == 0 ? nothing : n / count(support)
                                            for n in curves[k]],
                 "endpoint_identity_count" => count((labels[end] .== reference) .& support))
            for (k, (kind, label, support)) in enumerate(groups)]
    return (; transitions, returns, hidden, persistent,
            endpoint_identity = mask .& (labels[end] .== reference),
            sample_indices = indices, rows)
end

"""Binary-indicator sampled occupancy, including label 0 on the interior.

Outside-domain values are NaN, distinct from true zero occupancy.
"""
function species_occupancy(species, mask; stride::Int = 1)
    indices = 1:stride:length(species)
    out = fill(NaN, size(mask)..., 8)
    for s in 0:7
        counts = zeros(Int32, size(mask))
        for i in indices
            counts .+= species[i] .== s
        end
        layer = selectdim(out, ndims(out), s + 1)
        layer[mask] .= counts[mask] ./ length(indices)
    end
    return out
end

end
