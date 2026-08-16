# Deterministic tests of the radiation response direction.
# Pins the implemented semantics recorded in docs/audit_biofilms_potts.md §3:
# field maximal on the axis, sign(ΔH_rad) = sign(β_ion), negative-β cells drawn
# inward / radiosensitive species outward.

using Statistics: mean

# ---------- (a) Field shape: max on axis, monotone radial decay ----------

let
    p = SR.CPMParams(N = 20, n_cells_per_species = 1)
    st = SR.init_state(p; seed = 1)
    N = p.N
    zc = N ÷ 2
    yc = N ÷ 2          # y = 10 ⇒ y - N/2 = 0: pure-x radial line

    slice = st.radiation[:, :, zc]
    imax = Tuple(argmax(slice))
    @test SR.radial_dist(imax[1], imax[2], N) < 1.0   # global max on the axis

    radial_line = [st.radiation[x, yc, zc] for x in (N ÷ 2):N]
    @test issorted(radial_line, rev = true)            # monotone decay outward
    @test all(>(0), radial_line)                       # and dimensionless-positive
end

# ---------- (b) Exact local ΔH sign per species ----------
# Isolate the radiation term: J ≡ 0, λ_V = 0, melanin ≡ 0. Then for a copy that
# grows a single-site cell into a medium neighbor, ΔH = β_ion[sp] · I(target).

let
    p = SR.CPMParams(N = 20, n_cells_per_species = 1, λ_V = 0.0)
    st = SR.init_state(p; seed = 1)
    st.J .= 0
    fill!(st.melanin, 0.0)
    for i in eachindex(st.lattice)                     # wipe randomly-placed cells
        st.lattice[i] > 0 && (st.lattice[i] = Int32(0))
    end
    empty!(st.cells)

    # Source sites at two radii (both interior, R = 10): r = 0 (high I) and r = 7 (low I).
    hi_s, hi_t = (10, 10, 10), (11, 10, 10)
    lo_s, lo_t = (17, 10, 10), (18, 10, 10)
    I_hi = st.radiation[hi_t...]
    I_lo = st.radiation[lo_t...]
    @test I_hi > I_lo

    for sp in 1:SR.N_SPECIES
        st.cells[1] = SR.CellInfo(sp, 2, Float64[10, 10, 10])
        st.lattice[hi_s...] = Int32(1)
        st.lattice[lo_s...] = Int32(1)

        ΔH_hi = SR.compute_delta_H(st, hi_s..., hi_t...)
        ΔH_lo = SR.compute_delta_H(st, lo_s..., lo_t...)

        β = p.β_ion[sp]
        @test sign(ΔH_hi) == sign(β)                   # gain in high flux follows β
        @test sign(ΔH_lo) == sign(β)
        @test abs(ΔH_hi) > abs(ΔH_lo)                  # response decays with radius
        @test ΔH_hi ≈ β * I_hi atol = 1e-12            # nothing else contributes

        st.lattice[hi_s...] = Int32(0)                 # reset for next species
        st.lattice[lo_s...] = Int32(0)
        empty!(st.cells)
    end
end

# ---------- (c) Multi-seed drift with amplified synthetic β ----------
# Production β values are far too small to produce a stable short-run signal, so
# this uses a deliberately amplified SYNTHETIC coefficient (negative-signed for
# CN, radiosensitive-signed for SO, zero elsewhere) and averages across seeds.
# This tests the direction of the coupling, not calibrated biology.

mean_radius(st, sp, N) = begin
    rs = Float64[]
    for z in 1:N, y in 1:N, x in 1:N
        σ = st.lattice[x, y, z]
        if σ > 0 && SR.species_of(σ, st.cells) == sp
            push!(rs, SR.radial_dist(x, y, N))
        end
    end
    mean(rs)
end

let
    β_syn = zeros(7)
    β_syn[SR.CN] = -20.0        # amplified negative (toward-source) sign
    β_syn[SR.SO] = +20.0        # amplified radiosensitive sign
    N = 24
    n_mcs = 60
    seeds = 1:5

    drift_cn = Float64[]
    drift_so = Float64[]
    for seed in seeds
        p = SR.CPMParams(N = N, n_cells_per_species = 1, β_ion = β_syn)
        st = SR.init_state(p; seed)
        rng = SR.MersenneTwister(seed)
        r0_cn = mean_radius(st, SR.CN, N)
        r0_so = mean_radius(st, SR.SO, N)
        for _ in 1:n_mcs
            SR.mcs_step!(st, rng)
        end
        push!(drift_cn, mean_radius(st, SR.CN, N) - r0_cn)
        push!(drift_so, mean_radius(st, SR.SO, N) - r0_so)
    end

    @test mean(drift_cn) < 0    # negative-β cell moves toward the axis: a TROPISM
    @test mean(drift_so) > 0    # radiosensitive-signed cell moves outward
end
