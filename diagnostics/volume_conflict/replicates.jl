# Replicates PER CONFIGURATION.
#
# The earlier claim -- "s2.V > s1.V in 9 of 9, p = 0.002" -- observed nine (seed, ordering)
# configurations ONCE EACH and treated that as nine draws. On a backend that is not
# reproducible on identical input, one observation per configuration cannot separate a
# per-configuration effect from run-to-run noise. This runs R replicates of each of the nine.
#
# usage: julia --project=. diagnostics/volume_conflict/replicates.jl [R]
using Printf, Random, Statistics, LinearAlgebra
using JACC

const REPO = normpath(joinpath(@__DIR__, "..", ".."))
function load_serial()
    src = read(joinpath(REPO, "biofilms_potts.jl"), String)
    src = split(src, "#  13. Figure export")[1]
    M = Module(:SerialRef)
    Base.eval(M, :(using LinearAlgebra, Statistics, Random, Printf))
    Base.include_string(M, src, "biofilms_potts.jl")
    return M
end
const SR = load_serial()
try
    include(joinpath(REPO, "tests", "jacc_parity_tests.jl"))
catch
end

sd(v) = length(v) < 2 ? 0.0 : std(v)

function main(R)
    @printf("backend = %s   replicates per configuration = %d   V_MAX = %.3f\n\n",
            string(JACC.backend), R, V_MAX)
    @printf("%-16s %-17s %-17s %-17s %6s %6s %5s\n",
            "seed/ordering", "rate  mean±sd", "s1.V  mean±sd", "s2.V  mean±sd",
            "s1over", "s2over", "s2>s1")
    tot_s1over = 0; tot_s2over = 0; tot_up = 0; tot = 0
    per_config_up = Int[]
    for seed in (42, 43, 44), (name, ord) in PERMS
        rs = Float64[]; v1 = Float64[]; v2 = Float64[]
        for _ in 1:R
            A, E, _ = run_tables(seed, ord)
            s = pooled(A, E, 1:size(A, 1))
            h = size(A, 1) ÷ 2
            push!(rs, s.rate)
            push!(v1, pooled(A, E, 1:h).V)
            push!(v2, pooled(A, E, h+1:size(A, 1)).V)
        end
        s1o = count(>=(V_MAX), v1); s2o = count(>=(V_MAX), v2)
        up  = count(i -> v2[i] > v1[i], 1:R)
        tot_s1over += s1o; tot_s2over += s2o; tot_up += up; tot += R
        push!(per_config_up, up)
        @printf("%-16s %.5f±%.5f  %.5f±%.5f  %.5f±%.5f  %3d/%d %3d/%d %3d/%d\n",
                "$seed/$name", mean(rs), sd(rs), mean(v1), sd(v1), mean(v2), sd(v2),
                s1o, R, s2o, R, up, R)
    end
    println()
    @printf("  total runs                 : %d\n", tot)
    @printf("  s1.V >= V_MAX              : %d of %d (%.1f%%)\n", tot_s1over, tot, 100tot_s1over/tot)
    @printf("  s2.V >= V_MAX              : %d of %d (%.1f%%)\n", tot_s2over, tot, 100tot_s2over/tot)
    @printf("  s2.V > s1.V                : %d of %d (%.1f%%)\n", tot_up, tot, 100tot_up/tot)
    # exact two-sided binomial sign test against p = 0.5, no dependencies
    function binom_two_sided(k, n)
        pk(i) = Float64(binomial(big(n), big(i))) / 2.0^n
        target = pk(k) * (1 + 1e-12)
        sum(pk(i) for i in 0:n if pk(i) <= target)
    end
    @printf("  exact two-sided sign test  : p = %.4g\n", binom_two_sided(tot_up, tot))
    @printf("  configurations where s2>s1 in ALL %d replicates : %d of 9\n",
            R, count(==(R), per_config_up))
    @printf("  configurations where s2>s1 in NO replicate      : %d of 9\n",
            count(==(0), per_config_up))
end

main(length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 5)
