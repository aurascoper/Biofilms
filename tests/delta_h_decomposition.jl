# The four Hamiltonian terms, separately — and the proof that splitting them out
# changed nothing.
#
# `compute_delta_H` always computed ΔH_adh, ΔH_vol, ΔH_rad and ΔH_mel as four
# locals and discarded them on its return line. `compute_delta_H_terms` is that
# body with the sum removed; `compute_delta_H` now adds its fields back in the
# same order. The claim is bit-identity, and the claim is not free: summation
# order is OBSERVABLE here, so "same operands, same order" is doing work rather
# than describing an associativity that would hold anyway.
#
# contract_csv.jl covers the trajectory, and it has a floor: perturbing the
# melanin coefficient by one part in 5e6 leaves its CSV byte-identical, because
# a ~4e-8 shift in exp(-ΔH/T) has to straddle one of ~6.4M uniform draws to show
# up at all. That is fine for what it guards and useless for this: a refactor of
# the acceptance arithmetic needs an exact check, not a probabilistic one.

using Random   # runtests.jl brings in Test only; the sandbox module has its
                # own Random and this file must not borrow from it.

function _stepped_state(; seed = 7, mcs = 3)
    p = SR.CPMParams(N = 20, n_cells_per_species = 3)
    st = SR.init_state(p; seed = seed)
    rng = Random.Xoshiro(seed)
    # Melanin starts at zero everywhere and only grows through update_melanin!.
    # Without these steps ΔH_mel is 0.0 on every sampled move and the melanin
    # branch is dead surface a passing test would walk straight over.
    for _ in 1:mcs
        SR.mcs_step!(st, rng)
        SR.update_melanin!(st)
    end
    return p, st, rng
end

function _sample_moves(p, st, rng, n_draws)
    out = NTuple{6, Int}[]
    for _ in 1:n_draws
        sx, sy, sz = rand(rng, 1:p.N), rand(rng, 1:p.N), rand(rng, 1:p.N)
        st.interior[sx, sy, sz] || continue
        d = SR.NEIGHBORS_26[rand(rng, 1:26)]
        tx, ty, tz = sx + d[1], sy + d[2], sz + d[3]
        (1 <= tx <= p.N && 1 <= ty <= p.N && 1 <= tz <= p.N) || continue
        st.interior[tx, ty, tz] || continue
        st.lattice[sx, sy, sz] == st.lattice[tx, ty, tz] && continue
        push!(out, (sx, sy, sz, tx, ty, tz))
    end
    return out
end

@testset "ΔH decomposition" begin
    p, st, rng = _stepped_state()
    moves = _sample_moves(p, st, rng, 200_000)
    @test length(moves) > 1_000

    exact = 0
    reordered_differs = 0
    nonzero = zeros(Int, 4)
    for (sx, sy, sz, tx, ty, tz) in moves
        t = SR.compute_delta_H_terms(st, sx, sy, sz, tx, ty, tz)
        total = SR.compute_delta_H(st, sx, sy, sz, tx, ty, tz)
        # `===`, not `==`: -0.0 == 0.0 and NaN != NaN, and neither answer is
        # "the bits are the same".
        (t.adh + t.vol + t.rad + t.mel) === total && (exact += 1)
        (t.mel + t.rad + t.vol + t.adh) === total || (reordered_differs += 1)
        for (i, v) in enumerate((t.adh, t.vol, t.rad, t.mel))
            v != 0.0 && (nonzero[i] += 1)
        end
    end

    @test exact == length(moves)

    # THE CONTROL. If float addition were associative over these magnitudes the
    # bit-identity claim above would hold no matter what order the split summed
    # in, and this file would be asserting a tautology. It is not: adhesion and
    # volume are O(1), the radiation term is O(1e-5), and the melanin term is
    # O(0.5), so the small ones vanish or survive depending on where they land.
    @test reordered_differs > 0

    # SCOPE: the four branches of compute_delta_H_terms on the fixture below, not
# every path through the stepper.
# ALL FOUR BRANCHES MUST BE EXERCISED. A decomposition test that never sees
    # a nonzero melanin term proves nothing about the melanin term, and the
    # first version of this file did exactly that.
    for (i, name) in enumerate(("adh", "vol", "rad", "mel"))
        @test nonzero[i] > 0
    end
end

@testset "decisive labelling" begin
    T = 5.0
    mk(a, v, r, m) = (adh = a, vol = v, rad = r, mel = m)

    # --- the ΔH <= 0 branch: no draw was taken, u is NaN ---------------------
    # Removing the -2.0 adhesion term would put ΔH at +1.0, so the certainty is
    # gone -- but nothing here shows the move would have been REJECTED, and
    # saying so would need a number nobody drew.
    @test SR.decisive_label(mk(-2.0, 1.0, 0.0, 0.0), -1.0, NaN, T) ==
          SR.DRIVER_CONTINGENT
    # Every single removal still leaves ΔH <= 0: the move was never close.
    @test SR.decisive_label(mk(-5.0, -3.0, 0.0, 0.0), -8.0, NaN, T) ==
          SR.DRIVER_NONE

    # --- the ΔH > 0 branch: a draw exists, so nothing is contingent ----------
    # ΔH = 0.5, T = 5 -> threshold exp(-0.1) = 0.9048. u = 0.90 clears it.
    # Removing the -1.0 melanin term gives ΔH' = 1.5 -> exp(-0.3) = 0.7408,
    # which u does not clear. Melanin decided this move.
    @test SR.decisive_label(mk(1.5, 0.0, 0.0, -1.0), 0.5, 0.90, T) == 0x05
    # A term that HURT cannot be decisive: removing it only lowers ΔH.
    @test SR.decisive_label(mk(0.5, 0.0, 0.0, 0.0), 0.5, 0.90, T) ==
          SR.DRIVER_NONE
    # Two terms each individually decisive.
    @test SR.decisive_label(mk(2.5, 0.0, -1.0, -1.0), 0.5, 0.90, T) ==
          SR.DRIVER_MULTIPLE

    # --- THE DISJOINTNESS IS AN ASSERTION, NOT A REMARK ---------------------
    # If `contingent` ever appeared in the drawn branch, or a named term in the
    # undrawn one, the label would be answering two different questions under
    # one name and the colour key would be meaningless.
    rng = Random.Xoshiro(3)
    for _ in 1:20_000
        t = mk(randn(rng) * 2, randn(rng) * 2, randn(rng) * 1e-4,
               randn(rng) * 0.7)
        ΔH = t.adh + t.vol + t.rad + t.mel
        if ΔH <= 0
            @test SR.decisive_label(t, ΔH, NaN, T) in
                  (SR.DRIVER_NONE, SR.DRIVER_CONTINGENT)
        else
            u = rand(rng) * exp(-ΔH / T)      # accepted by construction
            lab = SR.decisive_label(t, ΔH, u, T)
            @test lab != SR.DRIVER_CONTINGENT
        end
    end

    # --- 0 IS NOT `none` ----------------------------------------------------
    d = SR.DriverCounts(4)
    SR.record_driver!(d, 2, 2, 2, SR.DRIVER_NONE)
    m = SR.modal_driver(d)
    @test m[2, 2, 2] == SR.DRIVER_NONE
    @test m[1, 1, 1] == 0x00          # never touched
    @test SR.n_accepted(d)[2, 2, 2] == 1
    @test SR.n_accepted(d)[1, 1, 1] == 0
end

@testset "the driver tally is inert" begin
    # THE INSTRUMENTATION MUST NOT MOVE THE TRAJECTORY. `u` is now kept rather
    # than discarded and the ternary is evaluated where an if/else was, so the
    # generator must still be consulted on exactly the moves it was before.
    # contract_csv.jl covers this at 100 MCS with its own floor; this covers it
    # directly, on the lattice itself, with and without the tally attached.
    function run(driver)
        p = SR.CPMParams(N = 16, n_cells_per_species = 3)
        st = SR.init_state(p; seed = 5)
        rng = Random.Xoshiro(5)
        for _ in 1:6
            SR.mcs_step!(st, rng; driver = driver === nothing ? nothing :
                                           SR.DriverCounts(p.N))
            SR.update_melanin!(st)
        end
        return copy(st.lattice), copy(st.melanin)
    end
    bare_lat, bare_mel = run(nothing)
    tally_lat, tally_mel = run(:on)
    @test bare_lat == tally_lat
    @test bare_mel == tally_mel

    # THE CONTROL: the comparison must be capable of failing. A different seed
    # must produce a different lattice, or the equality above is vacuous.
    p = SR.CPMParams(N = 16, n_cells_per_species = 3)
    other = SR.init_state(p; seed = 6)
    rng = Random.Xoshiro(6)
    for _ in 1:6
        SR.mcs_step!(other, rng)
        SR.update_melanin!(other)
    end
    @test other.lattice != bare_lat
end
