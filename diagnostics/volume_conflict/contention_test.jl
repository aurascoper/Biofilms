# Does one colour class actually contain proposals sharing a parcel id, and does the
# read timing change the outcome? Deterministic, host-only, no GPU.
using Test, Printf
include(joinpath(@__DIR__, "oracle.jl"))

const N   = 20
const SEED = UInt64(42)
const λV  = 10.0f0
const Vt  = Int32(120)
const T   = 5.0f0      # production default is T_cpm = 5.0f0; an earlier revision used 10.0f0

function fresh()
    lat, spec, vols, rad, mel, nut = init_host(N, 2, 42, 1.0, 2.0, 1.0)
    for idx in eachindex(mel); mel[idx] = Float32(0.01 * (idx % 7)); end
    return lat, spec, vols, rad, mel
end

J = build_J_matrix()

println("="^78)
println("1. Does a single colour class contain proposals that share a parcel id?")
println("="^78)
lat, spec, vols, rad, mel = fresh()

function survey_conflicts(lat, spec, vols, rad, mel, J)
  total_conflicts = 0
  for c in 0:7
    ox, oy, oz = c & 1, (c >> 1) & 1, (c >> 2) & 1
    l2, v2 = copy(lat), copy(vols)
    props = colour_pass!(l2, v2, spec, J, BETA_ION, MEL_COEF, rad, mel,
                         ox, oy, oz, SEED, UInt64(1), N, λV, Vt, T)
    con = contention(props)
    kinds = Dict{Symbol,Int}()
    for (_, _, _, kind) in con; kinds[kind] = get(kinds, kind, 0) + 1; end
    total_conflicts += length(con)
    @printf("  colour %d: %3d proposals, %4d conflicting pairs  %s\n",
            c, length(props), length(con), isempty(kinds) ? "" : string(kinds))
  end
  return total_conflicts
end
total_conflicts = survey_conflicts(lat, spec, vols, rad, mel, J)
println()
@testset "the checkerboard does not make parcel ids disjoint" begin
    @test total_conflicts > 0
end

println()
println("="^78)
println("2. Live vols vs sweep-start snapshot: same proposals, different decisions?")
println("="^78)
function compare_arms(lat, spec, vols, rad, mel, J)
  diverged = 0
  for c in 0:7
    ox, oy, oz = c & 1, (c >> 1) & 1, (c >> 2) & 1
    lA, vA = copy(lat), copy(vols)
    lB, vB = copy(lat), copy(vols)
    pA = colour_pass!(lA, vA, spec, J, BETA_ION, MEL_COEF, rad, mel,
                      ox, oy, oz, SEED, UInt64(1), N, λV, Vt, T; snapshot=false)
    pB = colour_pass!(lB, vB, spec, J, BETA_ION, MEL_COEF, rad, mel,
                      ox, oy, oz, SEED, UInt64(1), N, λV, Vt, T; snapshot=true)
    @assert length(pA) == length(pB) "arms proposed different sites"
    @assert all(i -> (pA[i].tx, pA[i].ty, pA[i].tz, pA[i].donor, pA[i].recipient) ==
                     (pB[i].tx, pB[i].ty, pB[i].tz, pB[i].donor, pB[i].recipient),
                1:length(pA)) "arms proposed different moves"
    d = count(i -> pA[i].accepted != pB[i].accepted, 1:length(pA))
    dh = count(i -> pA[i].dH != pB[i].dH, 1:length(pA))
    diverged += d
    @printf("  colour %d: %3d proposals, %3d differing ΔH, %3d differing decisions, accepted %3d vs %3d\n",
            c, length(pA), dh, d, count(p -> p.accepted, pA), count(p -> p.accepted, pB))
  end
  return diverged
end
diverged = compare_arms(lat, spec, vols, rad, mel, J)
println()
println("  total decisions changed by read timing: ", diverged)

println()
println("="^78)
println("3. Independent bookkeeping: does vols match the lattice histogram after a pass?")
println("="^78)
lC, vC = copy(lat), copy(vols)
function bookkeeping(lC, vC, spec, rad, mel, J)
  for c in 0:7
    ox, oy, oz = c & 1, (c >> 1) & 1, (c >> 2) & 1
    colour_pass!(lC, vC, spec, J, BETA_ION, MEL_COEF, rad, mel,
                 ox, oy, oz, SEED, UInt64(1), N, λV, Vt, T)
    h = id_histogram(lC, length(vC))
    ok = all(i -> h[i] == vC[i], eachindex(vC))
    @printf("  after colour %d: vols == lattice histogram : %s\n", c, ok ? "yes" : "NO")
  end
end
bookkeeping(lC, vC, spec, rad, mel, J)

@testset "sequential replay keeps vols consistent with the lattice" begin
    @test id_histogram(lC, length(vC)) == vC
end
