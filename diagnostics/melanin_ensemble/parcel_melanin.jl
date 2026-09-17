#!/usr/bin/env julia
# The two melanin estimators, computed from one visit over the lattice.
#
# `mel_site` gives every occupied site equal weight. A parcel of 200 sites counts 200
# times a parcel of 1, so the estimator is volume-weighted. `take_snapshot` in
# `biofilms_potts.jl` already reports exactly this, under the name `mean_melanin`.
#
# `mel_parcel` gives every parcel equal weight. It takes a per-cell mean first, then an
# unweighted mean over that species' cells. F. Fink's Odin port reports it. Before
# 2026-09-17 this repository computed it nowhere.
#
# THE TWO ESTIMATE DIFFERENT QUANTITIES. They agree only when every parcel of a species
# has the same site count. A control built from equal parcels would therefore prove
# nothing, so the controls in test_statistics.jl use unequal ones.
#
# THE DEFINITION IS THE PORT AUTHOR'S, QUOTED RATHER THAN PARAPHRASED (2026-09-16). A
# paraphrase of an observable is how the wrong quantity gets compared.
#   - Per occupied site with a live cell: cell_sum[slot] += melanin[i], cell_cnt[slot] += 1.
#   - Per live cell with count > 0: cell_mean = sum / cnt.
#   - Per species: unweighted mean of its live cells' cell_mean, 0 if none.
#     Dead and empty slots are skipped.
#
# This file takes arrays and a map rather than a CPMState. Its controls can then build a
# lattice by hand, and they need no model and no data file.

"""
    parcel_means(lattice, melanin, species_of, n_species)

Both melanin estimators per species, from one pass over `lattice`.

`lattice` holds 0 for medium, a negative value for wall, and a positive cell id
otherwise. `species_of` maps a live cell id to its species index. An id missing from that
map is skipped, which is the rule `take_snapshot` applies through `haskey(state.cells, s)`:
an id left on the lattice after its cell was removed belongs to neither mean.

Returns a named tuple of four vectors, indexed by species: `mel_site`, `mel_parcel`,
`n_sites` and `n_parcels`. A species with no occupied site gets 0.0, matching the model.

The visit order is Julia's own index order, which for an `Array{T,3}` moves x fastest and
z slowest. `take_snapshot` writes `for z in 1:N, y in 1:N, x in 1:N`, and that is the same
order, so the two site sums add their terms identically. The caller can therefore require
exact equality rather than a tolerance.
"""
function parcel_means(lattice, melanin, species_of::AbstractDict{Int,Int}, n_species::Int)
    size(lattice) == size(melanin) ||
        error("lattice $(size(lattice)) and melanin $(size(melanin)) have different shapes")
    n_species > 0 || error("n_species must be positive, got $n_species")

    cell_sum = Dict{Int,Float64}()
    cell_cnt = Dict{Int,Int}()
    site_sum = zeros(Float64, n_species)
    site_cnt = zeros(Int, n_species)

    @inbounds for i in eachindex(lattice)
        σ = Int(lattice[i])
        σ > 0 || continue
        haskey(species_of, σ) || continue
        s = species_of[σ]
        1 <= s <= n_species || error("cell $σ declares species $s, outside 1:$n_species")
        m = Float64(melanin[i])
        site_sum[s] += m
        site_cnt[s] += 1
        cell_sum[σ] = get(cell_sum, σ, 0.0) + m
        cell_cnt[σ] = get(cell_cnt, σ, 0) + 1
    end

    # A registered cell with no site on the lattice takes no part. The author's rule says
    # "per live cell with count > 0", and a zero count would divide by zero here.
    parcel_sum = zeros(Float64, n_species)
    parcel_n = zeros(Int, n_species)
    for (σ, cnt) in cell_cnt
        cnt > 0 || continue
        s = species_of[σ]
        parcel_sum[s] += cell_sum[σ] / cnt
        parcel_n[s] += 1
    end

    mel_site = [site_cnt[s] > 0 ? site_sum[s] / site_cnt[s] : 0.0 for s in 1:n_species]
    mel_parcel = [parcel_n[s] > 0 ? parcel_sum[s] / parcel_n[s] : 0.0 for s in 1:n_species]
    return (; mel_site, mel_parcel, n_sites = site_cnt, n_parcels = parcel_n)
end
