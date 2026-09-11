# Same input, same backend, twice. Does the acceptance table reproduce?
#
# This has to come before any statement attributing a statistic to a seed, an ordering or
# a half-window. A non-reproducible backend makes every such attribution a property of one
# sample rather than of the configuration.
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
    println("  (shipped assertions not executed here -- this driver reports)")
end

function repeat_same_input(n = 3; seed = 42, ord = IDENTITY)
    @printf("backend = %s   seed = %d   ordering = identity   replicates = %d\n\n",
            string(JACC.backend), seed, n)
    tabs = []
    for r in 1:n
        A, E, nf = run_tables(seed, ord)
        s = pooled(A, E, 1:size(A, 1))
        h = size(A, 1) ÷ 2
        s1 = pooled(A, E, 1:h); s2 = pooled(A, E, h+1:size(A, 1))
        push!(tabs, (A, E))
        @printf("  replicate %d: rate=%.5f  V=%.5f  s1.V=%.5f  s2.V=%.5f  %s%s\n",
                r, s.rate, s.V, s1.V, s2.V,
                s1.V >= V_MAX ? "FIRST-over " : "", s2.V >= V_MAX ? "SECOND-over" : "")
    end
    println()
    ident = all(i -> tabs[i][1] == tabs[1][1] && tabs[i][2] == tabs[1][2], 2:n)
    @printf("  all %d replicates byte-identical: %s\n", n, ident ? "YES" : "NO")
    if !ident
        d = sum(abs.(tabs[2][1] .- tabs[1][1]))
        @printf("  accepted-count L1 difference, replicate 1 vs 2: %d\n", d)
        @printf("  -> the backend is NOT reproducible on identical input, so any\n")
        @printf("     per-seed / per-half attribution describes one sample, not the run.\n")
    end
    return ident
end

repeat_same_input(3)
