#!/usr/bin/env julia
# Data-free tests for the sweep's statistics and its CSV contract.
#
#   julia diagnostics/melanin_ensemble/test_statistics.jl
#
# No committed sweep is required, no argument, no environment variable, so these assertions
# cannot be skipped by a missing sweep. The sign-test values below are exact rationals over
# 2^16 and are asserted as such rather than to a tolerance. The producer contract is tested
# by running sweep.jl itself for a few MCS and reading what it wrote through read_sweep.
using Test, Statistics
include(joinpath(@__DIR__, "analyse.jl"))
include(joinpath(@__DIR__, "figure_text.jl"))
include(joinpath(@__DIR__, "parcel_melanin.jl"))

const SWEEP = joinpath(@__DIR__, "sweep.jl")
sweep(out, args...) = success(pipeline(`$(Base.julia_cmd()) $SWEEP $out $args`; stdout = devnull, stderr = devnull))

@testset "melanin ensemble statistics" begin

@testset "exact two-sided binomial sign test" begin
    # Every outcome is no more likely than the median one, so the test returns 1.
    @test binom_two_sided(8, 16) ≈ 1.0
    # Only the two extremes are as unlikely as an extreme: 2/2^16.
    @test binom_two_sided(0, 16) ≈ 2 / 2^16
    @test binom_two_sided(16, 16) ≈ 2 / 2^16
    # 15 of 16: 2*(C(16,15) + C(16,16)) / 2^16 = 34/65536.
    @test binom_two_sided(15, 16) ≈ 34 / 65536
    # 12 of 16: 2*(1820 + 560 + 120 + 16 + 1) / 65536 = 5034/65536, which is ABOVE 0.05.
    @test binom_two_sided(12, 16) ≈ 5034 / 65536
    @test binom_two_sided(12, 16) > 0.05
    @test binom_two_sided(15, 16) < 0.05
    # Symmetric under k -> n-k, because the null is a fair coin.
    for k in 0:16
        @test binom_two_sided(k, 16) ≈ binom_two_sided(16 - k, 16)
    end
    # Monotone away from the centre.
    @test binom_two_sided(9, 16) > binom_two_sided(10, 16) > binom_two_sided(11, 16)
    # No observations: nothing is unlikely, p = 1.
    @test binom_two_sided(0, 0) ≈ 1.0
    # Past n = 1023 a Float64 2^n is Inf; a balanced result must still be p = 1 and an
    # extreme one small but nonzero.
    @test binom_two_sided(512, 1024) ≈ 1.0 atol = 1e-12
    @test 0 < binom_two_sided(0, 1024) < 1e-300
    # From n = 1076 one direction is below Float64's smallest subnormal; the value is kept
    # exact and printed from BigFloat, so it prints as a number and not as 0.
    @test binom_two_sided(1076, 1076) == 2 // big(2)^1076
    @test Float64(binom_two_sided(1076, 1076)) == 0.0            # what the old form printed
    @test fmtp(binom_two_sided(1076, 1076)) == "2.47e-324"
    @test fmtp(binom_two_sided(16, 16)) == "3.052e-05"           # the committed receipts' value
    # The pinned small-n values, to 1e-12, so the exact form cannot drift from the old one.
    @test isapprox(binom_two_sided(15, 16), 34 / 65536; atol = 1e-12)
    @test isapprox(binom_two_sided(12, 16), 5034 / 65536; atol = 1e-12)
    @test isapprox(binom_two_sided(0, 16), 2 / 65536; atol = 1e-12)
    @test isapprox(binom_two_sided(8, 16), 1.0; atol = 1e-12)
end

@testset "paired and unpaired spreads differ, and the direction is informative" begin
    # Two series that move together. The paired difference is constant, so its spread is
    # zero, while an independence assumption predicts sqrt(2) times the marginal sd. This
    # is the case where a pooled statistic overstates the uncertainty.
    a = [1.0, 2.0, 3.0, 4.0, 5.0]
    b = a .- 0.5
    d = a .- b
    @test std(d) ≈ 0.0 atol = 1e-12
    @test sqrt(var(a) + var(b)) > 2.2          # sqrt(2) * 1.5811
    @test mean(d) ≈ 0.5

    # Two series that move oppositely. The paired difference swings twice as far as
    # either series, so the paired spread EXCEEDS the independence prediction. This is
    # the case the melanin sweep is in, and the one where a pooled statistic understates.
    c = [5.0, 4.0, 3.0, 2.0, 1.0]
    dd = a .- c
    @test std(dd) > sqrt(var(a) + var(c))
    @test std(dd) ≈ 2 * std(a)
end

@testset "paired_stats, on synthetic seeds with a known answer" begin
    # Perfectly correlated: every seed shifts both species together, so the per-seed
    # difference is constant and the paired spread is zero. Independence would predict
    # sqrt(2) times the marginal sd, so a pooled statistic OVERSTATES here.
    rows = Dict(i => Dict("A" => Float64(i), "B" => Float64(i) - 0.5) for i in 1:5)
    r = paired_stats(rows, collect(1:5), "A", "B")
    @test r.mean ≈ 0.5
    @test r.sd_paired ≈ 0.0 atol = 1e-12
    @test r.sd_indep > r.sd_paired
    @test r.k == 5 && r.n == 5 && r.ties == 0
    @test r.p ≈ 2 / 2^5                       # 5 of 5 one way: 2/32

    # Anticorrelated: the difference swings twice as far as either series, so the paired
    # spread EXCEEDS the independence prediction. A pooled statistic UNDERSTATES here,
    # and this is the regime the melanin sweep is actually in.
    rows2 = Dict(i => Dict("A" => Float64(i), "B" => 6.0 - i) for i in 1:5)
    r2 = paired_stats(rows2, collect(1:5), "A", "B")
    @test r2.sd_paired > r2.sd_indep
    @test r2.sd_paired ≈ 2 * r2.sd_indep / sqrt(2)

    # A pair that is a coin flip gets a p-value that cannot reject anything.
    rows3 = Dict(1 => Dict("A"=>1.0,"B"=>0.0), 2 => Dict("A"=>0.0,"B"=>1.0),
                 3 => Dict("A"=>1.0,"B"=>0.0), 4 => Dict("A"=>0.0,"B"=>1.0))
    r3 = paired_stats(rows3, collect(1:4), "A", "B")
    @test r3.k == 2 && r3.p ≈ 1.0
    @test r3.mean ≈ 0.0 atol = 1e-12

    # The correlation is what drives the paired/independent gap, and its sign says which
    # way a pooled statistic errs. Perfectly correlated -> +1; anticorrelated -> -1.
    @test r.r ≈ 1.0 atol = 1e-12
    @test r2.r ≈ -1.0 atol = 1e-12
    @test r.sd_paired < r.sd_indep      # positive r: pooling overstates
    @test r2.sd_paired > r2.sd_indep    # negative r: pooling understates
    # var(hi - lo) = var(hi) + var(lo) - 2 cov, the identity the report rests on.
    let a2 = [1.0,2,3,4,5], b2 = [0.5,1.5,2.5,3.5,4.5]
        @test var(a2 .- b2) ≈ var(a2) + var(b2) - 2*cov(a2, b2) atol = 1e-12
    end

    # The separation is the mean over the POOLED sd, which is the statistic the paired
    # one replaces. Asserted so a change to either cannot silently swap them.
    @test r.separation ≈ r.mean / r.pooled
    @test r.pooled ≈ sqrt((var([1.0,2,3,4,5]) + var([0.5,1.5,2.5,3.5,4.5])) / 2)
end

@testset "a tie carries no sign and leaves the sign test" begin
    # All tied, the committed sweeps at MCS 0: the old test reported 0 of 16 with
    # p = 2/2^16, a directional result from data with no direction. Now: no trials, p = 1.
    tied = Dict(i => Dict("A" => 0.0, "B" => 0.0) for i in 1:16)
    t = paired_stats(tied, collect(1:16), "A", "B")
    @test t.k == 0 && t.n == 0 && t.ties == 16
    @test t.p ≈ 1.0
    # Three ties among five: the test sees the other two, both in one direction.
    mixed = Dict(1 => Dict("A"=>1.0,"B"=>1.0), 2 => Dict("A"=>2.0,"B"=>2.0),
                 3 => Dict("A"=>3.0,"B"=>3.0), 4 => Dict("A"=>5.0,"B"=>4.0),
                 5 => Dict("A"=>6.0,"B"=>4.0))
    m = paired_stats(mixed, collect(1:5), "A", "B")
    @test m.k == 2 && m.n == 2 && m.ties == 3
    @test m.p ≈ binom_two_sided(2, 2)
    # The mean and the spreads still see every seed; only the sign test drops ties.
    @test m.mean ≈ 3 / 5
    # Zero variance leaves the correlation and the separation undefined, not directional.
    @test isnan(t.r) && isnan(t.separation)
end

@testset "the sweep CSV contract, through the real reader" begin
    dir = mktempdir()
    header = "seed,mcs,species,species_name,alpha_M,volume,n_cells,mean_melanin"
    body = ["42,100,3,C. sphaerospermum,0.1400,230,2,1.5",
            "42,100,1,C. neoformans,0.1000,230,2,1.0",
            "42,100,5,A. niger,0.0650,230,2,0.5",
            "42,200,3,C. sphaerospermum,0.1400,230,2,9.9",   # another MCS, must be ignored
            "43,100,3,C. sphaerospermum,0.1400,230,2,1.6",
            "43,100,1,C. neoformans,0.1000,230,2,1.1",
            "43,100,5,A. niger,0.0650,230,2,0.4",
            "43,100,2,D. radiodurans,0.0000,230,2,7.7"]    # not a producer, must be ignored
    csv(name, head, rows) = (p = joinpath(dir, name); write(p, join(["# comment"; head; rows], "\n") * "\n"); p)

    rows, alpha = read_sweep(csv("ok.csv", header, body), 100)
    @test sort(collect(keys(rows))) == [42, 43]
    @test rows[42] == Dict("CS" => 1.5, "CN" => 1.0, "AN" => 0.5)
    @test rows[43]["AN"] == 0.4
    @test alpha == Dict("CS" => 0.14, "CN" => 0.1, "AN" => 0.065)

    # Columns are found by name, so a producer that reorders them still parses correctly.
    shuffled = "mean_melanin,alpha_M,species,mcs,seed"
    rows2, alpha2 = read_sweep(csv("shuffled.csv", shuffled, ["1.5,0.1400,3,100,42", "1.0,0.1000,1,100,42", "0.5,0.0650,5,100,42"]), 100)
    @test rows2[42] == Dict("CS" => 1.5, "CN" => 1.0, "AN" => 0.5)
    @test alpha2["CN"] == 0.1

    # A missing or renamed column refuses rather than parsing a neighbour as the observable.
    @test_throws ErrorException read_sweep(csv("renamed.csv", replace(header, "mean_melanin" => "melanin"), body), 100)
    @test_throws ErrorException read_sweep(csv("noalpha.csv", replace(header, "alpha_M" => "aM"), body), 100)
    @test_throws ErrorException read_sweep(csv("noheader.csv", "# only comments", String[]), 100)
    # A seed missing one producer used to pass the seed count and KeyError inside the pairing.
    @test_throws ErrorException read_sweep(csv("missing.csv", header, body[1:6]), 100)   # seed 43 lacks AN
    # A duplicated row is refused rather than silently overwriting the first.
    @test_throws ErrorException read_sweep(csv("dup.csv", header, [body; body[1]]), 100)
    # Coefficients in another order are refused by name: "the alpha_M direction" is CS > CN > AN.
    rev = [replace(replace(replace(l, ",0.1400," => ",X,"), ",0.0650," => ",0.1400,"), ",X," => ",0.0650,") for l in body]
    @test_throws ErrorException read_sweep(csv("reversed.csv", header, rev), 100)
    # A file whose early seeds carry other coefficients is refused at the first disagreement.
    @test_throws ErrorException read_sweep(csv("mixed.csv", header, [rev[1:3]; body[5:8]]), 100)

    # The three producers, by the species indices the model assigns.
    @test [id for (id, _, _) in PRODUCERS] == [3, 1, 5]
end

@testset "sweep.jl writes what read_sweep reads" begin
    # The producer itself, for four MCS on two seeds at the default N = 20: about three
    # seconds. This is the only test that touches biofilms_potts.jl, and it is the one that
    # binds the writer's header to the reader's expectations.
    dir = mktempdir()
    out = joinpath(dir, "tiny.csv")
    @test sweep(out, "--seeds", "42,43", "--mcs", "4", "--at", "4")
    rows, alpha = read_sweep(out, 4)
    @test sort(collect(keys(rows))) == [42, 43]
    @test all(haskey(rows[s], c) for s in (42, 43), c in ("CS", "CN", "AN"))
    @test alpha == Dict("CS" => 0.14, "CN" => 0.1, "AN" => 0.065)   # Table 2 midpoints
    @test all(v >= 0 for r in values(rows) for v in values(r))
    # The MCS 0 snapshot is in the file and is all ties; the reader must not confuse it.
    rows0, _ = read_sweep(out, 0)
    @test all(rows0[s]["CS"] == 0.0 for s in (42, 43))

    # Refusals happen before anything is written: a bad request leaves no file behind.
    @test !sweep(joinpath(dir, "late.csv"), "--seeds", "42", "--mcs", "4", "--at", "500")
    @test !isfile(joinpath(dir, "late.csv"))
    @test !sweep(joinpath(dir, "empty.csv"), "--seeds", "57:42", "--mcs", "4", "--at", "4")
    @test !isfile(joinpath(dir, "empty.csv"))
    @test !sweep(out, "--seeds", "42", "--mcs", "4", "--at", "4")       # destination exists
    @test !sweep(joinpath(dir, "dup.csv"), "--seeds", "42,42", "--mcs", "4", "--at", "4")
    @test !isfile(joinpath(dir, "dup.csv"))
end

@testset "figure_text: the figure reader refuses duplicates, and its strings follow the data" begin
    dir = mktempdir()
    meta = ["# melanin ensemble sweep", "# N=40 parcels_per_species=6 n_mcs=400 snapshot_interval=100"]
    header = "seed,mcs,species,species_name,alpha_M,volume,n_cells,mean_melanin"
    alpha_of(sp) = sp == 3 ? "0.1400" : sp == 1 ? "0.1000" : sp == 5 ? "0.0650" : "0.0000"
    row(seed, sp, mel) = "$seed,100,$sp,X,$(alpha_of(sp)),600,6,$mel"
    csv(name, rows) = (p = joinpath(dir, name); write(p, join([meta; header; rows], "\n") * "\n"); p)
    good = [row(42, 3, 1.4), row(42, 1, 0.9), row(42, 5, 0.8),
            row(43, 3, 1.3), row(43, 1, 1.0), row(43, 5, 0.7)]
    rows, alpha, m = read_at(csv("ok.csv", good), 100)
    @test sort(collect(keys(rows))) == [42, 43]
    @test rows[42][3] == 1.4 && alpha[3] == 0.14 && m["N"] == 40
    # A repeated (seed, species) measurement at this MCS is refused in either order --
    # equal or conflicting. Equal coefficients do not make two copies one run; the old
    # assignment kept whichever came last and said nothing.
    @test_throws ErrorException read_at(csv("dup_equal_last.csv", [good; good[1]]), 100)
    @test_throws ErrorException read_at(csv("dup_equal_first.csv", [good[1]; good]), 100)
    @test_throws ErrorException read_at(csv("dup_conflict_last.csv", [good; row(42, 3, 1.5)]), 100)
    @test_throws ErrorException read_at(csv("dup_conflict_first.csv", [row(42, 3, 1.5); good]), 100)
    # A repeat at another MCS is another observation, not a duplicate of this one.
    other, _, _ = read_at(csv("other_mcs.csv", [good; replace(good[1], ",100," => ",200,")]), 100)
    @test other[42][3] == 1.4
    # The mid-file alpha_M guard is kept beside the duplicate guard, not replaced by it.
    @test_throws ErrorException read_at(csv("alpha.csv", [good; replace(row(57, 3, 1.0), "0.1400" => "0.2000");
                                                           row(57, 1, 0.5); row(57, 5, 0.1)]), 100)

    # The title claims "every seed" only when every seed does.
    @test ordering_title(16, 16) == "Every seed descends: 16 of 16 display the α_M ordering"
    @test ordering_title(15, 16) == "15 of 16 display the α_M ordering"
    @test !occursin("Every seed", ordering_title(0, 2))
    # The provenance line prints the seeds that were plotted, not the interval between them.
    @test seed_set_label(42:57) == "42:57"
    @test seed_set_label([57, 42]) == "42,57"
    @test seed_set_label([42, 43, 44, 45, 50, 57]) == "42:45,50,57"
    @test seed_set_label([42]) == "42"
    @test_throws ErrorException seed_set_label(Int[])

    # The inverted, sparse fixture through the reader and the count figure.jl uses: seed 57
    # has CN above CS, and the seeds are 42 and 57 with nothing between.
    inverted = [row(42, 3, 1.4), row(42, 1, 0.9), row(42, 5, 0.8),
                row(57, 3, 1.0), row(57, 1, 1.2), row(57, 5, 0.7)]
    r2, _, _ = read_at(csv("inverted.csv", inverted), 100)
    seeds = sort(collect(keys(r2)))
    ordered = count(s -> r2[s][3] > r2[s][1] > r2[s][5], seeds)
    @test ordering_title(ordered, length(seeds)) == "1 of 2 display the α_M ordering"
    @test seed_set_label(seeds) == "42,57"
    # The all-ordered contiguous control keeps the committed wording byte for byte.
    @test ordering_title(2, 2) == "Every seed descends: 2 of 2 display the α_M ordering"
    @test seed_set_label([42, 43]) == "42:43"
    # Panel B: rank one is "the smallest", not "the smallest smallest".
    @test rank_title(1, 2) == "The published seed is the smallest of 2"
    @test rank_title(3, 16) == "The published seed is the third smallest of 16"
    @test rank_title(12, 16) == "The published seed is the 12th smallest of 16"
    @test rank_title(2, 16) == "The published seed is the second smallest of 16"
end

@testset "parcel_melanin: two estimators, and the case where they differ" begin
    # Two parcels of one species, of sizes 3 and 1. Site-weighted is (1+1+1+5)/4 = 2.0.
    # Parcel-weighted is (1.0 + 5.0)/2 = 3.0. That gap is the reason the estimator exists.
    lat = reshape(Int32[1, 1, 1, 2], 4, 1, 1)
    mel = reshape([1.0, 1.0, 1.0, 5.0], 4, 1, 1)
    r = parcel_means(lat, mel, Dict(1 => 1, 2 => 1), 1)
    @test r.mel_site[1] ≈ 2.0
    @test r.mel_parcel[1] ≈ 3.0
    @test r.n_sites[1] == 4
    @test r.n_parcels[1] == 2

    # EQUAL PARCELS MAKE THE TWO AGREE, so a control built from them would pass against
    # either estimator and prove nothing. The case above uses sizes 3 and 1 for that reason.
    eq = parcel_means(reshape(Int32[1, 1, 2, 2], 4, 1, 1),
                      reshape([1.0, 3.0, 5.0, 7.0], 4, 1, 1), Dict(1 => 1, 2 => 1), 1)
    @test eq.mel_site[1] ≈ 4.0
    @test eq.mel_parcel[1] ≈ 4.0

    # An id on the lattice with no registry entry belongs to neither mean. That is
    # take_snapshot's own rule, which tests haskey(state.cells, s) before it accumulates.
    stale = parcel_means(reshape(Int32[1, 1, 1, 9], 4, 1, 1), mel, Dict(1 => 1), 1)
    @test stale.mel_site[1] ≈ 1.0
    @test stale.mel_parcel[1] ≈ 1.0
    @test stale.n_sites[1] == 3

    # Medium is 0 and wall is negative, and neither takes part. A species with no site
    # reports 0.0, as the model does, rather than a NaN from a division by zero.
    none = parcel_means(reshape(Int32[0, -1, 2, 2], 4, 1, 1), mel, Dict(2 => 2), 2)
    @test none.mel_site[1] == 0.0
    @test none.mel_parcel[1] == 0.0
    @test none.n_parcels[1] == 0
    @test none.mel_site[2] ≈ 3.0        # sites 3 and 4 hold melanin 1.0 and 5.0

    # THIS EXPECTATION WAS WRONG THE FIRST TIME AND THE CODE WAS RIGHT. The draft asserted
    # 4.0 for the line above. Two sites holding 1.0 and 5.0 average to 3.0, and the
    # estimator said so. Recorded because a control that is merely green teaches nothing
    # about which side was checked.

    # Refusals: a melanin field of another shape, a species index outside the range, and a
    # species count of zero. Each would otherwise average the wrong set of sites.
    @test_throws ErrorException parcel_means(reshape(Int32[1], 1, 1, 1), mel, Dict(1 => 1), 1)
    @test_throws ErrorException parcel_means(lat, mel, Dict(1 => 3, 2 => 1), 1)
    @test_throws ErrorException parcel_means(lat, mel, Dict(1 => 1, 2 => 1), 0)
end

@testset "parcel_sweep.jl: the producer, its refusals, and the comparison's control" begin
    # RUN AS A SUBPROCESS, NOT INCLUDED. parcel_sweep.jl defines `COLUMNS` and `main`, and
    # analyse.jl, already included above this line, defines both names too. Including the
    # two would redefine a constant. The sweep testset runs its own producer this way.
    dir = mktempdir()
    PARCEL = joinpath(@__DIR__, "parcel_sweep.jl")
    parcel(args...) = success(pipeline(`$(Base.julia_cmd()) $PARCEL $args`;
                                       stdout = devnull, stderr = devnull))

    # N = 20 with two parcels is sweep.jl's own default size. Smaller than that the MODEL
    # refuses, not this file: at N = 8 with one parcel per species Random raises
    # "collection must be non-empty" inside mcs_step!, and sweep.jl fails there identically.
    out = joinpath(dir, "tiny.csv")
    @test parcel(out, "--seeds", "42,43", "--n", "20", "--parcels", "2", "--at", "4")

    lines = [l for l in eachline(out) if !startswith(l, "#") && !isempty(strip(l))]
    @test lines[1] == "seed,mcs,species,volume,ncells,mel_site,mel_parcel"
    @test length(lines) == 1 + 2 * 7          # the header, then seven species per seed
    body = [split(l, ",") for l in lines[2:end]]
    @test all(parse(Int, r[2]) == 4 for r in body)
    @test sort(unique(parse(Int, r[1]) for r in body)) == [42, 43]
    @test sort(unique(parse(Int, r[3]) for r in body)) == collect(1:7)
    @test all(parse(Float64, r[6]) >= 0 for r in body)
    @test all(parse(Float64, r[7]) >= 0 for r in body)
    # Six decimals, which is the precision the other implementation prints. A string
    # comparison between the two files means nothing unless both sides print the same way.
    @test all(length(split(String(r[7]), ".")[2]) == 6 for r in body)

    # Refusals, each before the model loads and before the destination opens.
    @test !parcel(out, "--seeds", "42", "--n", "20", "--parcels", "2", "--at", "4")
    for (name, args) in (("dup", ("--seeds", "42,42")),
                         ("empty", ("--seeds", "57:42")),
                         ("at0", ("--seeds", "42", "--at", "0")))
        p = joinpath(dir, "$name.csv")
        @test !parcel(p, args..., "--n", "20", "--parcels", "2")
        @test !isfile(p)
    end

    # THE COMPARISON'S OWN NEGATIVE CONTROL. A tool that cannot report a difference reads
    # exactly like one that always agrees, and the cross-implementation result rests on
    # this tool alone. One digit in the last decimal place has to be enough to refuse.
    @test parcel("--compare", out, out)
    row = lines[2]
    f = split(row, ",")
    digit = f[7][end] == '0' ? '1' : '0'
    bumped = join([f[1:6]; f[7][1:end-1] * string(digit)], ",")
    @test bumped != row
    mutated = joinpath(dir, "mutated.csv")
    write(mutated, replace(read(out, String), row => bumped))
    @test !parcel("--compare", mutated, out)

    # A file missing a column is refused rather than read by position.
    short = joinpath(dir, "short.csv")
    write(short, "seed,mcs,species,volume,ncells,mel_site\n42,4,3,1,1,1.0\n")
    @test !parcel("--compare", short, out)
end

end
