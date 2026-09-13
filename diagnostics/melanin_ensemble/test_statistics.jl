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

end
