# ---------- Checkerboard parity: the measurement nothing else makes ---------
#
# The port is an 8-color (2x2x2) checkerboard sweep. The characteristic defect
# of a checkerboard decomposition is a parity-correlated bias in accepted moves,
# and BOTH existing guards are blind to it: `jacc_port_tests.jl` compares
# kernels on identical inputs, so it passes when both kernels carry the same
# artifact, and `fixtures/serial_seed42.csv` pins the SERIAL stream, which has
# no sublattices at all.
#
# WHY A DISTRIBUTION AND NOT A MAP. d404438 refused a spatial map for the
# decisive-label work -- at 400 MCS only 6157 of 64000 voxels received any
# accepted move, median 3, so a mode over three samples is noise. The same
# argument applies here and the conclusion is stronger: a parity bias is a
# GLOBAL count comparison, so the map was never the instrument for it.
#
# WHY THE TABLE IS CONDITIONED ON OPPORTUNITY. A uniform null over eight
# classes assumes equal opportunity per class, and the kernel's early returns
# (wall, same-sigma, medium-into-medium, out-of-bounds) are all
# geometry-dependent -- an aggregate at an arbitrary lattice position has no
# reason to spread them evenly over a 2x2x2 decomposition. So the null "equal
# ACCEPTS per class" is false by construction and would report the shape of the
# domain as a decomposition artifact. `st` carries the denominator: a 2x8
# accepted/rejected contingency table asks about the acceptance RATE per class,
# which is the question.
#
# WHY RANDOM PERMUTATIONS AND NOT REVERSAL. `for c in color_order` is
# sequential and `vols` accumulates across passes (delta_H reads vols at :167
# while cpm_color! mutates it at :208 -- the documented `ponytail:` staleness),
# so the FIRST pass evaluates against sweep-start volumes and the LAST against
# volumes moved by seven passes. That is a deterministic, parity-correlated
# acceptance difference with nothing to do with the decomposition, and `c`
# indexes both spatial class and sequence position. Reversal maps
# position(c) = 7-c, which separates a MONOTONIC position effect but is blind to
# any position effect symmetric about the midpoint -- one transformation, one
# invariance. Random permutations drop the assumption about the effect's shape.
#
# WHY EFFECT SIZE AND NOT CHI-SQUARE. chi2 power scales with n and n here is
# 1e5-1e7 evaluated proposals, so a fraction-of-a-percent deviation clears any
# fixed critical value; a chi2 threshold would test "is there any asymmetry at
# all", and there always is. Worse, the cells are not independent -- an accepted
# move changes the lattice for every later pass -- so the nominal null does not
# apply and over-dispersion is expected from autocorrelation alone. Cramer's V
# and the max per-class rate deviation are the thresholds; chi2 and n are
# reported beside them and never asserted alone.
#
# basis_gate_ack: this file records NO radiodialysis quantity -- it records CPM
# acceptance counts. In the port, `delta_H` reads lat/vols/spec/J/beta/melc/rad/
# mel and never `nut`, and `melanin_k!` never reads it either, so the gated
# basis cannot reach an accepted move. Same claim as validate_serial.jl's, and
# the 10x-uptake testset below is what holds it to account.

using Random

const JACC_PARITY_PATH = joinpath(REPO, "biofilms_potts_jacc.jl")

# Spatial class DERIVED from target coordinates, never stored: deriving it is
# what keeps the class spatial when `color_order` permutes the pass sequence.
parity_class(tx, ty, tz) = (tx - 1) % 2 + 2 * ((ty - 1) % 2) + 4 * ((tz - 1) % 2)

# 2x8 accepted/rejected contingency statistics, conditioned on opportunity.
function parity_stats(a::Vector{Int}, e::Vector{Int})
    n = sum(e)
    p = n == 0 ? 0.0 : sum(a) / n
    χ2 = 0.0
    for c in 1:8, (o, ex) in ((a[c], e[c] * p), (e[c] - a[c], e[c] * (1 - p)))
        ex > 0 && (χ2 += (o - ex)^2 / ex)
    end
    rates = [e[c] == 0 ? NaN : a[c] / e[c] for c in 1:8]
    maxdev = p == 0 ? Inf : maximum(abs.(rates .- p)) / p
    return (; χ2, n, V = n == 0 ? Inf : sqrt(χ2 / n), rate = p, maxdev, rates)
end

const ParityPort = Module(:JaccParityPort)
Base.eval(ParityPort, :(using Random, Statistics, Printf))
Base.include_string(ParityPort, read(JACC_PARITY_PATH, String),
                    "biofilms_potts_jacc.jl")

# Separate top-level statement: `include` advances world age between them, so
# everything below can call into the freshly defined module directly.
const P_RP    = Base.eval(ParityPort, :RadiolysisParams)
const P_RUN   = Base.eval(ParityPort, :run_coupled)
const P_INIT  = Base.eval(ParityPort, :init_host)
const P_JMAT  = Base.eval(ParityPort, :build_J_matrix)
const P_KERN  = Base.eval(ParityPort, :cpm_color!)
const P_ARR   = Base.eval(ParityPort, :(JACC.array))
const P_HOST  = Base.eval(ParityPort, :(JACC.to_host))
const P_PFOR  = Base.eval(ParityPort, :(JACC.parallel_for))
const P_BETA  = Base.eval(ParityPort, :BETA_ION)
const P_MELC  = Base.eval(ParityPort, :MEL_COEF)

const IDENTITY = collect(0:7)
const PERMS = ["identity" => IDENTITY,
               "perm1"    => shuffle(MersenneTwister(1), collect(0:7)),
               "perm2"    => shuffle(MersenneTwister(2), collect(0:7))]

# OBSERVED, NOT DERIVED. T_cpm = 5.0f0 gives no derivable band, and a bound read
# off the quantity's domain -- an acceptance rate is a probability, so
# `0 < rate < 1` cannot fail for the reason it exists -- would report green at
# 1e-9. These come from the spread over 3 seeds x 3 orderings at N=20/50 MCS on
# the threads backend (V 0.0050-0.0111, maxdev 0.017-0.049, rate 0.171-0.191),
# widened ~2x. They are settings-specific: at N=40 the rate is ~0.05, because
# V_target is fixed while lattice volume is not.
const V_MAX      = 0.025
const MAXDEV_MAX = 0.12
const RATE_BAND  = (0.14, 0.23)

mk_rp(; kw...) = P_RP(; Nr = 40, Ddot_R = 1.0, c_ext = 1.0,
                      basis_gate_ack = true, kw...)

# Per-sweep (accepted, evaluated) per spatial class, plus the finiteness guard
# on the discriminator.
function run_tables(seed, order; N = 20, n_mcs = 50, rp = mk_rp())
        A = zeros(Int, n_mcs, 8); E = zeros(Int, n_mcs, 8); nonfinite = 0
        cb = (mcs, st, dh) -> begin
            nonfinite += count(!isfinite, @view dh[st .!= 0])
            for tz in 1:N, ty in 1:N, tx in 1:N
                s = st[tx, ty, tz]; s == 0 && continue
                c = parity_class(tx, ty, tz) + 1
                E[mcs, c] += 1; s == UInt8(2) && (A[mcs, c] += 1)
            end
        end
    P_RUN(; seed, n_mcs, N, rp, verbose = false,
          color_order = order, on_sweep = cb)
    return A, E, nonfinite
end

pooled(A, E, rows) = parity_stats(vec(sum(A[rows, :], dims = 1)),
                                  vec(sum(E[rows, :], dims = 1)))

@testset "JACC checkerboard parity" begin
    @testset "the class encode/decode pair agrees with the kernel" begin
        # The kernel encodes c -> (ox,oy,oz) -> tx = 2*(i-1)+ox+1; this file
        # decodes coordinates back to a class. Two mappings maintained apart,
        # and the class cannot be read back from the kernel because deriving it
        # is what makes it survive permutation. So assert the round trip against
        # the REAL kernel: one color pass in isolation, every site it touched
        # must decode to exactly the color that was passed.
        N = 20; Nh = N ÷ 2
        lat_h, spec_h, vols_h, rad_h, mel_h, _ =
            P_INIT(N, 6, 42, 1.0, 2.0, 1.0)
        lat = P_ARR(lat_h); vols = P_ARR(vols_h); spec = P_ARR(spec_h)
        J = P_ARR(P_JMAT())
        βv = P_ARR(P_BETA); melc = P_ARR(P_MELC)
        rad = P_ARR(rad_h); mel = P_ARR(mel_h)
        st = P_ARR(zeros(UInt8, N, N, N)); dh = P_ARR(zeros(Float32, N, N, N))
        for c in 0:7
            fill!(st, UInt8(0))
            P_PFOR((Nh, Nh, Nh), P_KERN,
                lat, vols, spec, J, βv, melc, rad, mel,
                Int32(c & 1), Int32((c >> 1) & 1), Int32((c >> 2) & 1),
                UInt64(7), UInt64(c), Int32(N), 10.0f0, Int32(120), 5.0f0, st, dh)
            st_h = P_HOST(st)
            touched = [(tx, ty, tz) for tz in 1:N, ty in 1:N, tx in 1:N
                       if st_h[tx, ty, tz] != 0]
            @test !isempty(touched)
            @test all(t -> parity_class(t...) == c, touched)
        end
    end

    @testset "color_order must be a permutation of 0:7" begin
        # Without this guard [0,0,1,2,3,4,5,6] double-updates one class, never
        # updates another, and yields a silently wrong simulation.
        for bad in ([0, 0, 1, 2, 3, 4, 5, 6], collect(0:6), collect(1:8))
            @test_throws AssertionError P_RUN(; N = 20, n_mcs = 1, seed = 42,
                rp = mk_rp(), verbose = false, color_order = bad)
        end
    end

    @testset "per-sweep reset, discriminator, and all eight classes" begin
        A, E, nonfinite = run_tables(42, IDENTITY; n_mcs = 3)
        # A per-COLOR-PASS reset would leave only the last color populated and
        # would look perfectly clean; this is what catches it.
        @test all(>(0), sum(E, dims = 1))
        # `dh` carries no sentinel, so it must be finite wherever `st` says a
        # proposal was evaluated -- otherwise "never proposed" and "NaN" collapse
        # and every denominator below silently shrinks.
        @test nonfinite == 0
    end

    @testset "the guard can fire" begin
        # A threshold nobody has seen fail is not a guard. Synthetic table with
        # one class held 20% below the others, at the n this suite actually
        # reaches: it must clear BOTH thresholds.
        e = fill(17_000, 8); a = round.(Int, e .* 0.185); a[4] = round(Int, 17_000 * 0.148)
        s = parity_stats(a, e)
        @test s.V > V_MAX
        @test s.maxdev > MAXDEV_MAX
    end

    @testset "no decomposition artifact across seeds and orderings" begin
        results = Pair{String,Any}[]
        for seed in (42, 43, 44), (name, ord) in PERMS
            A, E, _ = run_tables(seed, ord)
            s = pooled(A, E, 1:size(A, 1))
            push!(results, "seed $seed $name" => s)

            @test RATE_BAND[1] < s.rate < RATE_BAND[2]
            @test s.V < V_MAX
            @test s.maxdev < MAXDEV_MAX

            # Stationarity before any pooled number is believed: an
            # initialization transient and a persistent bias are different
            # findings and only the second is a decomposition artifact.
            h = size(A, 1) ÷ 2
            s1 = pooled(A, E, 1:h); s2 = pooled(A, E, h+1:size(A, 1))
            @test s1.V < V_MAX && s2.V < V_MAX
        end

        # Reported, not asserted: below V_MAX the imbalance is not attributable
        # to either spatial class or pass position, and the permuted run is a
        # DIFFERENT trajectory rather than the same system observed differently,
        # so a third outcome -- matching neither -- is reachable and is what the
        # measured data shows. If a threshold above ever fails, these two rows
        # are what attributes it: a bias tracking spatial class under every
        # ordering is a real decomposition artifact; one tracking pass position
        # is the `vols` staleness and is not a finding about the checkerboard.
        for (label, s) in results
            @info "parity" label n=s.n rate=round(s.rate, digits=5) V=round(s.V, digits=5) χ2=round(s.χ2, digits=1) maxdev=round(s.maxdev, digits=4) rates=round.(s.rates, digits=4)
        end
    end

    @testset "the exemption's claim is true: acceptance does not see the basis" begin
        # basis_gate_ack above is acknowledged on the grounds that acceptance
        # counts cannot depend on the gated biomass basis. That is a CLAIM. The
        # basis enters only through `uptake`/`nut`, so 10x the uptake constants
        # must leave the contingency table byte-identical. The day acceptance
        # reads the nutrient field, this fails and the exemption is re-argued.
        base = run_tables(42, IDENTITY; n_mcs = 10)
        pert = run_tables(42, IDENTITY; n_mcs = 10,
                          rp = mk_rp(k_ads = 0.5, k_red = 0.2))
        @test base[1] == pert[1]
        @test base[2] == pert[2]
    end
end
