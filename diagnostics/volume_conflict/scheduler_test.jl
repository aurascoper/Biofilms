# The scheduler's only claim: identical to the sequential oracle, every proposal,
# every ΔH, every decision, and the final state -- over a full trajectory, not one pass.
using Test, Printf
include(joinpath(@__DIR__, "scheduler.jl"))

const N = 20; const SEED = UInt64(42)
const λV = 10.0f0; const Vt = Int32(120); const T = 5.0f0

function fresh()
    lat, spec, vols, rad, mel, _ = init_host(N, 2, 42, 1.0, 2.0, 1.0f0)
    for idx in eachindex(mel); mel[idx] = Float32(0.01 * (idx % 7)); end
    return lat, spec, vols, rad, mel
end

key(p) = (p.tx, p.ty, p.tz)

function run_both(n_mcs)
    J = build_J_matrix()
    latA, spec, volsA, rad, mel = fresh()
    latB, _,    volsB, _,   _   = fresh()
    mismatch = 0; total = 0
    batches = Int[]; sizes = Int[]
    for mcs in 1:n_mcs, c in 0:7
        ox, oy, oz = c & 1, (c >> 1) & 1, (c >> 2) & 1
        step = UInt64(mcs * 8 + c)
        pA = colour_pass!(latA, volsA, spec, J, BETA_ION, MEL_COEF, rad, mel,
                          ox, oy, oz, SEED, step, N, λV, Vt, T)
        pB, nb, nsites = scheduled_pass!(latB, volsB, spec, J, BETA_ION, MEL_COEF, rad, mel,
                                         ox, oy, oz, SEED, step, N, λV, Vt, T)
        push!(batches, nb); push!(sizes, nsites)
        dA = Dict(key(p) => p for p in pA)
        dB = Dict(key(p) => p for p in pB)
        total += length(pA)
        length(pA) == length(pB) || (mismatch += 1; continue)
        keys(dA) == keys(dB) || (mismatch += 1; continue)
        for k in keys(dA)
            a, b = dA[k], dB[k]
            (a.donor == b.donor && a.recipient == b.recipient &&
             a.dH == b.dH && a.draw == b.draw && a.accepted == b.accepted) || (mismatch += 1)
        end
    end
    return latA, volsA, latB, volsB, mismatch, total, batches, sizes
end

println("Scheduled batches vs the sequential oracle, 50 MCS x 8 colours\n")
latA, volsA, latB, volsB, mismatch, total, batches, sizes = run_both(50)

@printf("  proposals compared          : %d\n", total)
@printf("  proposals differing         : %d\n", mismatch)
@printf("  final lattice identical     : %s\n", latA == latB ? "yes" : "NO")
@printf("  final volumes identical     : %s\n", volsA == volsB ? "yes" : "NO")
println()
@printf("  batches per pass  min/mean/max : %d / %.1f / %d\n",
        minimum(batches), sum(batches)/length(batches), maximum(batches))
@printf("  proposals per pass mean        : %.1f\n", sum(sizes)/length(sizes))
@printf("  mean parallelism (props/batch) : %.2f\n", (sum(sizes)/length(sizes)) / (sum(batches)/length(batches)))

@testset "the scheduler reproduces the sequential oracle exactly" begin
    @test mismatch == 0
    @test latA == latB
    @test volsA == volsB
end

# A conflict-free partition that does NOT preserve order must NOT be assumed equivalent.
println("\nControl: is order preservation load-bearing, or is conflict-freedom enough?")
function shuffled_batches(n_mcs)
    J = build_J_matrix()
    latC, spec, volsC, rad, mel = fresh()
    latD, _,    volsD, _,   _   = fresh()
    diff = 0
    for mcs in 1:n_mcs, c in 0:7
        ox, oy, oz = c & 1, (c >> 1) & 1, (c >> 2) & 1
        step = UInt64(mcs * 8 + c)
        colour_pass!(latC, volsC, spec, J, BETA_ION, MEL_COEF, rad, mel,
                     ox, oy, oz, SEED, step, N, λV, Vt, T)
        # reverse the within-class enumeration order: still id-disjoint per batch when
        # re-layered, but conflicting pairs now resolve in the opposite order
        sites = reverse(propose(latD, ox, oy, oz, SEED, step, N))
        b = schedule(sites)
        frozen = copy(volsD); cur = 0
        for idx in sortperm(b)
            s = sites[idx]
            if b[idx] != cur; frozen = copy(volsD); cur = b[idx]; end
            ΔH = delta_H(latD, frozen, spec, J, BETA_ION, MEL_COEF, rad, mel,
                         s.tx, s.ty, s.tz, s.donor, s.recipient, N, λV, Vt)
            if ΔH <= 0.0f0 || u01(s.r2) < exp(-ΔH / T)
                s.recipient > Int32(0) && (volsD[s.recipient] -= Int32(1))
                s.donor > Int32(0) && (volsD[s.donor] += Int32(1))
                latD[s.tx, s.ty, s.tz] = s.donor
            end
        end
    end
    return latC == latD, volsC == volsD
end
sl, sv = shuffled_batches(50)
@printf("  reversed-order batching reproduces the oracle: lattice %s, volumes %s\n",
        sl ? "yes" : "NO", sv ? "yes" : "NO")
println("  (if NO, conflict-freedom alone is insufficient and order preservation is required)")
