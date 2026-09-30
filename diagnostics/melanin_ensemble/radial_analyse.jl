#!/usr/bin/env julia
# Paired analysis of a radial ensemble, against the tests registered before it was run.
#
#   julia diagnostics/melanin_ensemble/radial_analyse.jl <radial.csv> [--at 100]
#
# THE COMPARISON IS PAIRED WITHIN A SEED, across ARMS. Every arm of a seed starts from
# the same placement (radial_sweep.jl checks it on every seed), so the per-seed difference
# between two arms is the effect of the term that was switched off, on that placement.
# The statistics are the ones analyse.jl uses for species within a seed: the paired mean
# and spread, the within-run correlation, and the exact sign test with ties dropped.
#
# THE PRIMARY TESTS AND THEIR THRESHOLDS ARE CONSTANTS HERE AND PROSE IN
# PREDICTION_radial_n40_p6.md, and test_statistics.jl asserts the two agree. They were
# written before the ensemble ran; the file says so and says what was seen beforehand.
#
#   MOVES     under full - no_mel, for CN and CS: sign test p < ALPHA (two-sided) AND
#             |mean delta r/R| >= MEOI.
#   EQUIVALENT under full - no_rad, for every species: the 90% CI of the paired mean
#             (two one-sided tests at ALPHA) lies inside +-EQUIV_BOUND.
#
# Every other species x contrast cell is reported and labelled exploratory.

using Printf, Statistics
include(joinpath(@__DIR__, "stats.jl"))

const ALPHA       = 0.05
const MEOI        = 0.05     # r/R; one lattice unit at N = 40, below a parcel's own radius
const EQUIV_BOUND = 0.05     # r/R; the same unit, as the bound for "does not move"
const PRIMARY_MOVE_SPECIES = (1, 3)          # CN, CS: the two species with beta_ion < 0
const CONTRASTS = (("full", "no_mel"), ("full", "no_rad"), ("full", "null"), ("no_rad", "null"))

const RADIAL_COLUMNS = ("seed", "arm", "mcs", "species", "species_name", "mean_r_over_R", "mean_melanin")

"""
Read a radial_sweep.jl CSV at one MCS into species => seed => arm => r/R, plus the
species names, the header metadata, and the same structure for the MCS 0 seeding frame.

Every seed must carry every (arm, species) pair the file's header declares; a seed
missing one used to slip past a seed count and fail inside the pairing.
"""
function read_radial(path::AbstractString, at::Int)
    rr   = Dict{Int, Dict{Int, Dict{String, Float64}}}()
    rr0  = Dict{Int, Dict{Int, Dict{String, Float64}}}()
    names = Dict{Int, String}()
    meta = Dict{String, String}()
    col = nothing
    for line in eachline(path)
        isempty(strip(line)) && continue
        if startswith(line, "#")
            for m in eachmatch(r"(\w+)=(\S+)", line)
                meta[m[1]] = m[2]
            end
            continue
        end
        f = split(line, ",")
        if isnothing(col)
            idx = Dict(strip(String(f[i])) => i for i in eachindex(f))
            absent = [c for c in RADIAL_COLUMNS if !haskey(idx, c)]
            isempty(absent) || error("$path: header lacks $(join(absent, ", ")): $line")
            col = Dict(c => idx[c] for c in RADIAL_COLUMNS)
            continue
        end
        mcs = parse(Int, f[col["mcs"]])
        target = mcs == at ? rr : mcs == 0 ? rr0 : nothing
        isnothing(target) && continue
        sp   = parse(Int, f[col["species"]])
        seed = parse(Int, f[col["seed"]])
        arm  = String(strip(f[col["arm"]]))
        names[sp] = String(strip(f[col["species_name"]]))
        s = get!(get!(target, sp, Dict{Int, Dict{String, Float64}}()), seed, Dict{String, Float64}())
        haskey(s, arm) && error("$path: duplicate row for seed $seed arm $arm species $sp MCS $mcs")
        s[arm] = parse(Float64, f[col["mean_r_over_R"]])
    end
    isnothing(col) && error("$path: no header line")
    isempty(rr) && error("$path: no rows at MCS $at")
    haskey(meta, "arms") || error("$path: header does not declare arms=")
    arms = String.(split(meta["arms"], ","))
    for (label, d) in (("MCS $at", rr), ("MCS 0", rr0)), (sp, per_seed) in d, (seed, per_arm) in per_seed
        for arm in arms
            haskey(per_arm, arm) || error("$path: seed $seed species $sp has no $arm row at $label")
        end
    end
    # The pairing rests on a shared placement. The producer checks it exactly; the reader
    # checks it again on the printed values, because a file assembled from two runs would
    # pass the producer's check in each and fail it here.
    for (sp, per_seed) in rr0, (seed, per_arm) in per_seed
        vals = unique(per_arm[a] for a in arms)
        length(vals) == 1 || error("$path: seed $seed species $sp: the MCS 0 seeding frame differs " *
                                   "across arms ($(join(vals, ", "))); the arms do not share a placement")
    end
    return (at = rr, seeding = rr0, names = names, arms = arms, meta = meta)
end

"Paired statistics of arm `hi` minus arm `lo` for one species, over the seeds both carry."
function arm_contrast(per_seed, hi::AbstractString, lo::AbstractString)
    seeds = sort(collect(keys(per_seed)))
    paired_stats(per_seed, seeds, hi, lo)
end

"Two-sided normal quantile for the CI the equivalence test uses; 1.6449 at 90%."
const Z90 = 1.6448536269514722

"The 90% CI of a paired mean, and whether it lies inside +-bound (two one-sided tests)."
function equivalence(r, n::Int, bound::Float64)
    half = Z90 * r.sd_paired / sqrt(n)
    lo, hi = r.mean - half, r.mean + half
    (lo = lo, hi = hi, inside = lo > -bound && hi < bound)
end

"The pre-registered verdict for a MOVES cell."
moves(r) = r.p < ALPHA && abs(r.mean) >= MEOI

function radial_main(args)
    isempty(args) && error("usage: radial_analyse.jl <radial.csv> [--at 100]")
    path = args[1]
    i = findfirst(==("--at"), args)
    at = isnothing(i) ? 100 : parse(Int, args[i + 1])
    d = read_radial(path, at)
    species = sort(collect(keys(d.at)))
    seeds = sort(collect(keys(d.at[species[1]])))
    n = length(seeds)
    n >= 2 || error("paired statistics need at least two seeds; $n at MCS $at in $path")
    for c in CONTRASTS, a in c
        a in d.arms || error("$path declares arms $(join(d.arms, ",")); contrast $(c[1]) - $(c[2]) needs $a")
    end
    for sp in PRIMARY_MOVE_SPECIES
        sp in species || error("$path carries no rows for species $sp, which primary test P1 names")
    end

    @printf("radial ensemble %s, MCS %d, %d seeds, arms %s\n", basename(path), at, n, join(d.arms, ","))
    @printf("produced by julia %s on %s\n", get(d.meta, "julia_version", "?"), get(d.meta, "machine", "?"))
    @printf("thresholds: alpha %.2f two-sided, MEOI %.2f r/R, equivalence bound +-%.2f r/R (90%% CI)\n\n",
            ALPHA, MEOI, EQUIV_BOUND)

    println("seeding frame (MCS 0, identical across arms), mean r/R over seeds:")
    for sp in species
        x = [first(values(d.seeding[sp][s])) for s in seeds]
        @printf("  %-22s  mean %.4f   sd %.4f\n", d.names[sp], mean(x), std(x))
    end

    println("\nper arm, marginal mean r/R at MCS $at:")
    @printf("  %-22s", "")
    foreach(a -> @printf("  %14s", a), d.arms); println()
    for sp in species
        @printf("  %-22s", d.names[sp])
        for a in d.arms
            x = [d.at[sp][s][a] for s in seeds]
            @printf("  %6.4f +-%5.4f", mean(x), std(x))
        end
        println()
    end

    println("\nPRIMARY: full - no_mel, CN and CS (pre-registered: MOVES if p < alpha and |mean| >= MEOI)")
    for sp in PRIMARY_MOVE_SPECIES
        r = arm_contrast(d.at[sp], "full", "no_mel")
        @printf("  %-22s  mean %+.4f  paired sd %.4f  %d of %d seeds positive, p = %s  ->  %s\n",
                d.names[sp], r.mean, r.sd_paired, r.k, r.n, fmtp(r.p), moves(r) ? "MOVES" : "does not meet the criterion")
    end

    println("\nPRIMARY: full - no_rad, every species (pre-registered: EQUIVALENT if the 90% CI is inside the bound)")
    for sp in species
        r = arm_contrast(d.at[sp], "full", "no_rad")
        e = equivalence(r, n, EQUIV_BOUND)
        @printf("  %-22s  mean %+.4f  90%% CI [%+.4f, %+.4f]  ->  %s\n",
                d.names[sp], r.mean, e.lo, e.hi, e.inside ? "EQUIVALENT within the bound" : "NOT shown equivalent")
    end

    println("\nEXPLORATORY: every species x contrast, paired within seed")
    for (hi, lo) in CONTRASTS
        @printf("  %s - %s\n", hi, lo)
        for sp in species
            r = arm_contrast(d.at[sp], hi, lo)
            corr = isnan(r.r) ? "undefined" : @sprintf("%+.3f", r.r)
            @printf("    %-22s  mean %+.4f  paired sd %.4f  within-run r %s  %d of %d positive, p = %s%s\n",
                    d.names[sp], r.mean, r.sd_paired, corr, r.k, r.n, fmtp(r.p),
                    r.ties == 0 ? "" : @sprintf(" (%d ties)", r.ties))
        end
    end

    # The manuscript's one number: the C. neoformans - B. subtilis separation, seed 42,
    # was about 0.15 r/R. Here it is per arm over every seed, and its full - null difference.
    if 1 in species && 4 in species
        println("\nCN - BS separation within each arm (the single-seed figure's quantity), paired within seed:")
        sepd = Dict(s => Dict(a => d.at[1][s][a] - d.at[4][s][a] for a in d.arms) for s in seeds)
        for a in d.arms
            x = [sepd[s][a] for s in seeds]
            nz = filter(!=(0.0), x)
            @printf("  %-8s  mean %+.4f  sd %.4f  %d of %d seeds positive, p = %s\n",
                    a, mean(x), std(x), count(>(0), nz), length(nz), fmtp(binom_two_sided(count(>(0), nz), length(nz))))
        end
        r = paired_stats(sepd, seeds, "full", "null")
        @printf("  full - null of that separation: mean %+.4f  paired sd %.4f  %d of %d positive, p = %s\n",
                r.mean, r.sd_paired, r.k, r.n, fmtp(r.p))
    end

    println("\nThe arms are declared inputs. These numbers say whether switching a term off moves the")
    println("mean radius the model computes; they say nothing about organisms.")
end

# Parenthesised: `@__FILE__ && x` parses as `@__FILE__(&& x)`.
if abspath(PROGRAM_FILE) == @__FILE__
    radial_main(ARGS)
end
