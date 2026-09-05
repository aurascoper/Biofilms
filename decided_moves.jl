#!/usr/bin/env julia
# Table 4 (Share of accepted moves reversed by removing a single Hamiltonian
# term), reproduced from the shipped code.
#
#   julia --project=. decided_moves.jl            # Table 4, both columns
#   julia --project=. decided_moves.jl 42 43 44   # per-seed at 400 MCS
#
# WHY THIS FILE EXISTS. Table 4 was published from a configuration NO ENTRY
# POINT REPRODUCED. `run_simulation` and `run_simulation_coupled` both seed
# MersenneTwister and give 14281 accepted moves at 100 MCS against the table's
# 16037; the table came from the idiom in tests/delta_h_decomposition.jl --
# `mcs_step!` driven by hand with Random.Xoshiro(seed) and `update_melanin!`
# each sweep -- and neither b39fb8a nor d404438 committed a harness or recorded
# the parameters. A published table whose run cannot be re-executed is a number
# nobody can check, which is what PP-62-11 records and what this discharges.
#
# It uses only the shipped API: `mcs_step!` already takes the `driver` keyword
# and `DriverCounts` already accumulates the labels. Nothing here recomputes a
# label or reimplements the counterfactual -- `decisive_label` is called by
# `mcs_step!`, not by this file, so the table cannot drift from the rule that
# produces it.
#
# NO RADIODIALYSIS AND THEREFORE NO BASIS GATE. The CPM trajectory does not read
# the nutrient field (`compute_delta_H_terms` takes lattice, volumes, species,
# J, beta_ion, melanin and radiation, and never `nut`), which is why the coupled
# and uncoupled paths give identical accepted counts. This file steps the CPM
# and the melanin field only, so it acknowledges no gated basis.

using Printf, Random

const HERE = @__DIR__

# Same split-marker load as validate_serial.jl: the serial monolith minus its
# CairoMakie section, so this runs without a plotting stack.
function load_serial()
    src = read(joinpath(HERE, "biofilms_potts.jl"), String)
    src = split(src, "#  13. Figure export")[1]
    M = Module(:SerialRef)
    Base.eval(M, :(using LinearAlgebra, Statistics, Random, Printf))
    Base.include_string(M, src, "biofilms_potts.jl")
    return M
end

"""
    decided_moves(SR, seed, n_mcs; N, n_cells_per_species) -> Vector{Int}

Accepted moves by decisive-term label, in `DRIVER_LABELS` order.
"""
function decided_moves(SR::Module, seed::Int, n_mcs::Int;
                       N::Int = 40, n_cells_per_species::Int = 6)
    p   = SR.CPMParams(; N, n_cells_per_species)
    st  = SR.init_state(p; seed)
    rng = Random.Xoshiro(seed)          # NOT MersenneTwister -- see the header
    d   = SR.DriverCounts(p.N)
    for _ in 1:n_mcs
        SR.mcs_step!(st, rng; driver = d)
        SR.update_melanin!(st)
    end
    return vec(sum(d.counts, dims = (1, 2, 3)))
end

# The published table, so this file states what it must reproduce rather than
# printing numbers nobody can check against anything.
const PUBLISHED = Dict(
    (42, 100) => [6184, 1207,  710, 0,  2,  14,  7920],
    (42, 400) => [35372, 3296, 2127, 0, 23,  68, 12717],
    (43, 400) => [52605, 2455, 1995, 0, 20,  75, 11315],
    (44, 400) => [64439, 3041, 2044, 1, 43, 130, 14276],
)

function report(SR, seed, n_mcs)
    labels = SR.DRIVER_LABELS
    c   = decided_moves(SR, seed, n_mcs)
    acc = sum(c)
    ref = get(PUBLISHED, (seed, n_mcs), nothing)
    @printf("\nseed %d, %d MCS -- %d accepted moves\n", seed, n_mcs, acc)
    @printf("  %-12s %8s %8s   %s\n", "term removed", "count", "share", "published")
    for k in eachindex(labels)
        @printf("  %-12s %8d %7.2f%%   %s\n", labels[k], c[k], 100 * c[k] / acc,
                ref === nothing ? "-" :
                c[k] == ref[k] ? "ok" : "MISMATCH (was $(ref[k]))")
    end
    return ref === nothing || c == ref
end

function main(args)
    SR = load_serial()
    ok = true
    if isempty(args)
        ok &= Base.invokelatest(report, SR, 42, 100)   # Table 4, left column
        ok &= Base.invokelatest(report, SR, 42, 400)   # Table 4, right column
    else
        for s in parse.(Int, args)
            ok &= Base.invokelatest(report, SR, s, 400)
        end
    end
    println()
    if ok
        println("reproduces the published table")
    else
        println("DOES NOT reproduce the published table -- see MISMATCH above")
        exit(1)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
