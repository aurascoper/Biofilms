# Report s1.V and s2.V SEPARATELY for every seed x ordering.
#
# jacc_parity_tests.jl:209 asserts `s1.V < V_MAX && s2.V < V_MAX`. Julia prints no operand
# values for a `&&` compound, so a failure names neither the offending half nor its value --
# which is why the three half-window failures on Metal could not be attributed from the log.
#
# Reuses the shipped run_tables / pooled / parity_stats / PERMS so the statistics are the
# file's own, not a reimplementation. Asserts nothing; reports.
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

# The shipped testsets run on include and throw on a divergent backend. Everything this
# diagnostic needs is defined above them, so the definitions survive the throw.
try
    include(joinpath(REPO, "tests", "jacc_parity_tests.jl"))
catch e
    println("  include threw: ", first(split(string(e), "\n")))
    println("  (shipped assertions not executed here -- this driver reports, it does not assert)")
end

function halves()
    @printf("backend = %s   V_MAX = %.4f   MAXDEV_MAX = %.3f   RATE_BAND = (%.2f, %.2f)\n\n",
            string(JACC.backend), V_MAX, MAXDEV_MAX, RATE_BAND[1], RATE_BAND[2])
    @printf("%-18s %8s %8s %8s %10s %10s %8s\n",
            "seed / ordering", "rate", "V", "maxdev", "s1.V", "s2.V", "over")
    nr = 0; np = 0; n1 = 0; n2 = 0; n209 = 0
    for seed in (42, 43, 44), (name, ord) in PERMS
        A, E, _ = run_tables(seed, ord)
        s  = pooled(A, E, 1:size(A, 1))
        h  = size(A, 1) ÷ 2
        s1 = pooled(A, E, 1:h)
        s2 = pooled(A, E, h+1:size(A, 1))
        br = !(RATE_BAND[1] < s.rate < RATE_BAND[2])
        bp = (s.V >= V_MAX) || (s.maxdev >= MAXDEV_MAX)
        b1 = s1.V >= V_MAX
        b2 = s2.V >= V_MAX
        nr += br; np += bp; n1 += b1; n2 += b2; n209 += (b1 || b2)
        @printf("%-18s %8.5f %8.5f %8.4f %10s %10s %8s\n",
                "$seed / $name", s.rate, s.V, s.maxdev,
                (b1 ? "*" : " ") * @sprintf("%.5f", s1.V),
                (b2 ? "*" : " ") * @sprintf("%.5f", s2.V),
                b1 && b2 ? "both" : b1 ? "first" : b2 ? "second" : "-")
    end
    println()
    @printf("  rate outside band        : %d of 9\n", nr)
    @printf("  pooled V or maxdev over  : %d of 9\n", np)
    @printf("  FIRST  half V >= V_MAX   : %d of 9\n", n1)
    @printf("  SECOND half V >= V_MAX   : %d of 9\n", n2)
    @printf("  runs failing line 209    : %d of 9\n", n209)
end
halves()
