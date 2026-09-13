#!/usr/bin/env julia
# Paired analysis of a melanin seed sweep.
#
#   julia diagnostics/melanin_ensemble/analyse.jl <sweep.csv> [--at 100]
#
# The comparison is PAIRED WITHIN A SEED, and that is the point of this file. The three
# melanin producers share one lattice in every run, so a pooled standard deviation across
# species treats 16 runs as 32 independent draws and misstates the spread: it understates
# when the two series are negatively correlated within a run and overstates when they are
# positively correlated, and the committed sweeps show one pair each way. The paired
# statistic needs no independence assumption; the ordering count needs no distributional
# assumption at all. All three are reported side by side so the difference is visible
# rather than asserted.
#
# alpha_M is an input. Nothing here measures it, and agreement with it is a statement
# about how reliably one trajectory DISPLAYS an input, not about whether the input is right.
using Printf, Statistics

const PRODUCERS = [(3, "CS", "C. sphaerospermum"), (1, "CN", "C. neoformans"), (5, "AN", "A. niger")]

"Exact two-sided binomial sign test at p = 1/2: the total probability of every outcome no more likely than the one observed."
function binom_two_sided(k::Int, n::Int)
    pk(i) = Float64(binomial(big(n), big(i))) / 2.0^n
    target = pk(k) * (1 + 1e-12)
    sum(pk(i) for i in 0:n if pk(i) <= target)
end

"""
    paired_stats(rows, seeds, hi, lo)

Difference statistics for one species pair, computed WITHIN each seed.

`sd_paired` is the spread of the per-seed difference and assumes nothing. `sd_indep` is
what that spread would be if the two series were independent draws. When `sd_paired`
exceeds `sd_indep` the two are negatively correlated within a run, and a pooled statistic
understates the uncertainty; when it falls below, the pooled statistic overstates it.
Extracted from `main` so this, the load-bearing computation in this file, is reachable
from the data-free tests. A mutation check found it unreachable while it was inline.
"""
function paired_stats(rows, seeds, hi::AbstractString, lo::AbstractString)
    d = [rows[s][hi] - rows[s][lo] for s in seeds]
    xh = [rows[s][hi] for s in seeds]
    xl = [rows[s][lo] for s in seeds]
    pooled = sqrt((var(xh) + var(xl)) / 2)
    # The within-run correlation is what makes the paired and independent spreads differ:
    # var(hi - lo) = var(hi) + var(lo) - 2 cov(hi, lo). Reported because the sign of it
    # says which way a pooled statistic is wrong.
    r = cor(xh, xl)
    # A tie carries no sign, so it leaves the sign test entirely: out of k AND out of n.
    # Counting it in n alone reads every tie as a vote against, and an all-tie pair
    # (MCS 0, where nothing has been produced yet) came back as 0 of 16 with p = 3e-5
    # instead of the no-information answer, p = 1.
    nz = filter(!=(0.0), d)
    k, n = count(>(0), nz), length(nz)
    (mean = mean(d), sd_paired = std(d), sd_indep = sqrt(var(xh) + var(xl)),
     pooled = pooled, separation = mean(d) / pooled, r = r,
     k = k, n = n, ties = length(d) - n, p = binom_two_sided(k, n))
end

# The columns the analysis needs, found BY NAME in the producer's header. A reordering in
# sweep.jl is then harmless, and a missing or renamed column refuses instead of parsing
# the wrong field as the observable.
const COLUMNS = ("seed", "mcs", "species", "alpha_M", "mean_melanin")

function read_sweep(path, at)
    rows = Dict{Int, Dict{String, Float64}}()
    alpha = Dict{String, Float64}()
    col = nothing
    for line in eachline(path)
        (isempty(line) || startswith(line, "#")) && continue
        f = split(line, ",")
        if isnothing(col)
            idx = Dict(name => i for (i, name) in enumerate(f))
            absent = [c for c in COLUMNS if !haskey(idx, c)]
            isempty(absent) || error("$path: header lacks $(join(absent, ", ")): $line")
            col = Dict(c => idx[c] for c in COLUMNS)
            continue
        end
        parse(Int, f[col["mcs"]]) == at || continue
        sp = parse(Int, f[col["species"]])
        for (id, code, _) in PRODUCERS
            sp == id || continue
            seed = parse(Int, f[col["seed"]])
            get!(rows, seed, Dict{String, Float64}())[code] = parse(Float64, f[col["mean_melanin"]])
            alpha[code] = parse(Float64, f[col["alpha_M"]])
        end
    end
    isnothing(col) && error("$path: no header line")
    rows, alpha
end

function main(args)
    isempty(args) && error("usage: analyse.jl <sweep.csv> [--at 100]")
    path = args[1]
    i = findfirst(==("--at"), args)
    at = isnothing(i) ? 100 : parse(Int, args[i + 1])
    rows, alpha = read_sweep(path, at)
    seeds = sort(collect(keys(rows)))
    n = length(seeds)
    # One seed passes a non-empty check and hands std, var and cor a single observation,
    # which prints NaN spreads under a heading that says "paired".
    n >= 2 || error("paired statistics need at least two seeds; $n at MCS $at in $path")

    @printf("sweep %s, MCS %d, %d seeds\n\n", basename(path), at, n)
    @printf("%6s %9s %9s %9s   %s\n", "seed", "CS", "CN", "AN", "ordering")
    for s in seeds
        v = rows[s]
        tag = (v["CS"] > v["CN"] > v["AN"]) ? "alpha_M" :
              join(filter(!isempty, [v["AN"] > v["CN"] ? "AN>CN" : "",
                                     v["CN"] > v["CS"] ? "CN>CS" : ""]), " ")
        @printf("%6d %9.3f %9.3f %9.3f   %s\n", s, v["CS"], v["CN"], v["AN"], tag)
    end

    println("\nper species, marginal:")
    for (_, code, name) in PRODUCERS
        x = [rows[s][code] for s in seeds]
        @printf("  %-3s alpha_M %.3f   mean %.3f   sd %.3f   range %.3f - %.3f   %s\n",
                code, alpha[code], mean(x), std(x), minimum(x), maximum(x), name)
    end

    agree = count(s -> rows[s]["CS"] > rows[s]["CN"] > rows[s]["AN"], seeds)
    @printf("\nseeds displaying the full alpha_M ordering: %d of %d\n", agree, n)

    println("\npairwise, paired within seed:")
    for (hi, lo) in (("CS", "CN"), ("CN", "AN"))
        r = paired_stats(rows, seeds, hi, lo)
        @printf("  %s - %s   mean %+.3f\n", hi, lo, r.mean)
        @printf("      paired sd of the difference   %.3f   <- needs no independence assumption\n", r.sd_paired)
        @printf("      if the two were independent   %.3f   (pooled sd %.3f, separation %.2f)\n",
                r.sd_indep, r.pooled, r.separation)
        @printf("      within-run correlation        %+.3f   -> pooling %s the uncertainty\n",
                r.r, r.sd_paired > r.sd_indep ? "UNDERSTATES" : "overstates")
        @printf("      %d of %d seeds in the alpha_M direction, exact two-sided sign test p = %.4g%s\n",
                r.k, r.n, r.p, r.ties == 0 ? "" : @sprintf(" (%d ties dropped)", r.ties))
    end

    println("\nThe ordering is an input. These numbers describe how reliably one trajectory")
    println("displays it, and say nothing about whether the input is correct.")
end

# Parenthesised: `@__FILE__ && x` parses as `@__FILE__(&& x)`.
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
