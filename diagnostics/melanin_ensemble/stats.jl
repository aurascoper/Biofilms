# Statistics shared by analyse.jl and radial_analyse.jl: the exact sign test and the paired
# difference. Factored out of analyse.jl on 2026-09-24 so the radial ensemble could use the
# same three functions without a copy; the tests in test_statistics.jl reach them
# through either includer. Nothing here reads a file.
using Printf, Statistics

"""
Exact two-sided binomial sign test at p = 1/2: the total probability of every outcome no
more likely than the one observed, as an exact rational. Not converted to Float64: the
earlier form divided each term by 2.0^n, Inf past n = 1023, and the one after it converted
the exact result, which underflows to 0 from n = 1076 in one direction. `fmtp` prints it.
"""
function binom_two_sided(k::Int, n::Int)
    c = [binomial(big(n), big(i)) for i in 0:n]
    sum(c[i] for i in eachindex(c) if c[i] <= c[k + 1]) // big(2)^n
end

"Format a probability without losing it to Float64 underflow."
fmtp(p) = @sprintf("%.4g", BigFloat(p))

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
