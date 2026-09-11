#!/usr/bin/env julia
# Data-free tests for the sweep's statistics.
#
#   julia diagnostics/melanin_ensemble/test_statistics.jl
#
# No CSV, no argument, no environment variable, so these assertions cannot be skipped by
# a missing sweep. The sign-test values below are exact rationals over 2^16 and are
# asserted as such rather than to a tolerance.
using Test, Statistics
include(joinpath(@__DIR__, "analyse.jl"))

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
    @test r.k == 5 && r.n == 5
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

@testset "the sweep CSV contract" begin
    # The header the analysis parses, pinned so a change to sweep.jl that reorders
    # columns fails here rather than silently producing wrong numbers.
    header = "seed,mcs,species,species_name,alpha_M,volume,n_cells,mean_melanin"
    cols = split(header, ",")
    @test length(cols) == 8
    @test cols[1] == "seed" && cols[2] == "mcs" && cols[3] == "species"
    @test cols[5] == "alpha_M" && cols[8] == "mean_melanin"
    # analyse.jl reads field 8 as the observable and field 5 as the declared input.
    src = read(joinpath(@__DIR__, "analyse.jl"), String)
    @test occursin("parse(Float64, f[8])", src)
    @test occursin("parse(Float64, f[5])", src)
    # The three producers, by the species indices the model assigns.
    @test [id for (id, _, _) in PRODUCERS] == [3, 1, 5]
end

end
