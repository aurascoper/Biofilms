#!/usr/bin/env julia
# The per-proposal ΔH_rad statistics quoted in Section 6.2, reproduced from the
# shipped code.
#
#   julia --project=. rad_proposals.jl            # the three Table 4 runs
#   julia --project=. rad_proposals.jl 42         # one seed
#
# WHY THIS FILE EXISTS. Section 6.2 states that 26.7% of evaluated proposals
# carry ΔH_rad < -5e-5, over 3968838 proposals across seeds 42, 43 and 44 at
# 400 MCS. VERSION 1.2 STATED 29.8% OVER A RUN IT DID NOT NAME, and nothing in
# the repository could produce either number: `decided_moves.jl` counts accepted
# moves by decisive label and never touches evaluated proposals, and the
# measurement behind the restatement was an in-memory rewrite of
# `biofilms_potts.jl` that was never committed. A corrected number nobody can
# re-run is the defect PP-62-11 records, committed inside the fix for it.
# Raised by Codex on pull request #23.
#
# It uses the `on_proposal` hook on `mcs_step!`, which is called for every
# evaluated proposal before the acceptance draw and never consults the
# generator. tests/rad_proposals_tests.jl pins that the hook does not move the
# trajectory, with a different-seed control so the comparison can fail.
#
# NO RADIODIALYSIS AND THEREFORE NO BASIS GATE, for the reason decided_moves.jl
# gives: compute_delta_H_terms never reads the nutrient field.

using Printf, Random

const HERE = @__DIR__

function load_serial()
    src = read(joinpath(HERE, "biofilms_potts.jl"), String)
    src = split(src, "#  13. Figure export")[1]
    M = Module(:SerialRef)
    Base.eval(M, :(using LinearAlgebra, Statistics, Random, Printf))
    Base.include_string(M, src, "biofilms_potts.jl")
    return M
end

"The threshold §6.2 reports against: the bound version 1.1 claimed."
const THRESH = -5e-5

"""
    rad_proposals(SR, seed, n_mcs; N, n_cells_per_species) -> NamedTuple

`n` evaluated proposals, `below` of them carrying ΔH_rad < THRESH, and the
extremes of ΔH_rad over the run.

Counters rather than a stored vector: the three runs evaluate about four million
proposals and holding them costs more than the simulation.
"""
function rad_proposals(SR::Module, seed::Int, n_mcs::Int;
                       N::Int = 40, n_cells_per_species::Int = 6)
    p   = SR.CPMParams(; N, n_cells_per_species)
    st  = SR.init_state(p; seed)
    rng = Random.Xoshiro(seed)          # matches decided_moves.jl, NOT MersenneTwister
    n = Ref(0); below = Ref(0)
    lo = Ref(Inf); hi = Ref(-Inf)
    hook = function (terms, _ΔH)        # never touches `rng` -- see the header
        n[] += 1
        terms.rad < THRESH && (below[] += 1)
        lo[] = min(lo[], terms.rad)
        hi[] = max(hi[], terms.rad)
    end
    for _ in 1:n_mcs
        SR.mcs_step!(st, rng; on_proposal = hook)
        SR.update_melanin!(st)
    end
    return (; n = n[], below = below[], frac = below[] / n[], min = lo[], max = hi[])
end

# Measured 2026-08-29 at N=40, 6 parcels per species, 400 MCS, Xoshiro(seed).
# `min` is the acceptance-favouring extreme; `max` is its mirror, and at seeds 43
# and 44 it reaches +0.07505 -- the pairwise extremum, in the disfavouring
# direction. The acceptance-favouring -0.07505 is NOT reached at N=40, which is
# what §6.2 says and what makes the single-role bound survive its own run.
# THESE COUNTS ARE MEASURED. The first draft of this table back-computed them
# from the reported fractions -- 1287796 * 0.2783 -- which is a number that looks
# like evidence and is not one. Running this file printed the real counts and
# they disagreed in the third digit.
const PUBLISHED = Dict(
    42 => (n = 1_287_796, below =  358_404, min = -0.071342207),
    43 => (n = 1_340_171, below =  344_591, min = -0.061028003),
    44 => (n = 1_340_871, below =  357_192, min = -0.075),
)
const POOLED_N    = 3_968_838
const POOLED_FRAC = 0.267

function main(args)
    seeds = isempty(args) ? [42, 43, 44] : parse.(Int, args)
    SR = load_serial()
    tn = 0; tb = 0; lo = Inf; hi = -Inf
    @printf("%-6s %12s %12s %8s %14s %14s\n",
            "seed", "evaluated", "< -5e-5", "frac", "min dH_rad", "max dH_rad")
    for s in seeds
        # invokelatest for the same reason decided_moves.jl uses it: the module
        # is built at run time and its bindings do not exist in this world yet.
        r = Base.invokelatest(rad_proposals, SR, s, 400)
        @printf("%-6d %12d %12d %8.4f %14.8g %14.8g\n", s, r.n, r.below, r.frac, r.min, r.max)
        tn += r.n; tb += r.below; lo = min(lo, r.min); hi = max(hi, r.max)
    end
    @printf("%-6s %12d %12d %8.4f %14.8g %14.8g\n", "pooled", tn, tb, tb / tn, lo, hi)
    @printf("\nSection 6.2 states %.1f%% over %d evaluated proposals.\n",
            100 * POOLED_FRAC, POOLED_N)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
