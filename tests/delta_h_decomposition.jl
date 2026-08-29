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

    # ALL FOUR BRANCHES MUST BE EXERCISED. A decomposition test that never sees
    # a nonzero melanin term proves nothing about the melanin term, and the
    # first version of this file did exactly that.
    for (i, name) in enumerate(("adh", "vol", "rad", "mel"))
        @test nonzero[i] > 0
    end
end
