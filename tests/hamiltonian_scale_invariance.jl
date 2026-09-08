# The Metropolis rule reads only ΔH/T (biofilms_potts.jl:800-802):
#
#     u = ΔH <= 0 ? NaN : rand(rng)
#     accept = ΔH <= 0 ? true : u < exp(-ΔH / p.T_cpm)
#
# Each active term is linear in exactly one coefficient -- ΔH_adh in J, ΔH_vol in
# λ_V (:549), ΔH_rad in β_ion (:566) -- so scaling every coefficient AND T_cpm by a
# common c > 0 leaves every acceptance decision unchanged, including which branch of
# the ternary is taken, since `ΔH <= 0` is sign-preserving. The generator is therefore
# consulted at exactly the same proposals and the stream is untouched.
#
# CONSEQUENCE: the shipped T_cpm = 5.0 is one redundant degree of freedom. No fit
# against a trajectory can identify the scale of the coefficient vector. This is a
# third non-identifiability alongside the melanin product α_M·k (NEWS-MEL-03) and the
# pitch degeneracy V_physical = V_sites·a³ (docs/calibration/cpm_spatial_calibration.md).
#
# With c a power of two the claim is BIT-EXACT, not statistical: scaling by a power of
# two is exact in IEEE-754 and commutes with rounding through the sums and the
# division, so (2ΔH)/(2T) is exactly ΔH/T. That lets this be an equality rather than a
# tolerance -- the strongest form available.

using Random

const SI_N, SI_NC = 12, 2
const SI_P = SR.CPMParams(N = SI_N, n_cells_per_species = SI_NC)

# A second copy of the serial source with the hard-coded melanin literal DOUBLED in
# memory. The substitution is taken from the artifact and its application is asserted
# below: a mutation that silently fails to apply reports "no difference" and reads
# exactly like the result it was meant to test.
const SI_MEL_SUBS = Ref(0)
const SI_SR_MEL2 = let
    src = split(read(joinpath(REPO, "biofilms_potts.jl"), String), "#  13. Figure export")[1]
    patched = replace(src, "0.5 * M_local" => "1.0 * M_local")
    SI_MEL_SUBS[] = count("1.0 * M_local", patched) - count("1.0 * M_local", src)
    M = Module(:SerialMelaninDoubled)
    Base.eval(M, :(using LinearAlgebra, Statistics, Random, Printf))
    Base.include_string(M, patched, "biofilms_potts.jl")
    M
end

"A state under module `M` with every H coefficient scaled by cH and T by cT."
function si_state(M, cH, cT; melanin = 0.0, seed = 7)
    s = M.init_state(M.CPMParams(N = SI_N, n_cells_per_species = SI_NC); seed = seed)
    s.J .*= cH
    s.params = M.CPMParams(N = SI_N, n_cells_per_species = SI_NC,
                           T_cpm = cT * SI_P.T_cpm, λ_V = cH * SI_P.λ_V,
                           β_ion = cH .* SI_P.β_ion)
    s.melanin .= melanin
    s
end

function si_trace(cH, cT; melanin = 0.0, n = 30, seed = 7)
    s = si_state(SR, cH, cT; melanin = melanin, seed = seed)
    rng = MersenneTwister(4242)
    for _ in 1:n
        SR.mcs_step!(s, rng)
    end
    ids = sort(collect(keys(s.cells)))
    (s.lattice, ids, [s.cells[i].volume for i in ids])
end

si_same(a, b) = a[1] == b[1] && a[2] == b[2] && a[3] == b[3]

"""
The 1-based MCS at which two runs first differ, or 0 if they never do within n.

A bare "they diverged" is weak: two runs that differ at all will differ everywhere
eventually, so the count of differing sites says little. WHEN they first differ says
whether the cause bites immediately or is nearly cancelling.
"""
function si_first_divergence(MA, MB, cHa, cTa, cHb, cTb; melanin = 0.0, n = 30)
    a = si_state(MA, cHa, cTa; melanin = melanin)
    b = si_state(MB, cHb, cTb; melanin = melanin)
    ra, rb = MersenneTwister(4242), MersenneTwister(4242)
    for k in 1:n
        MA.mcs_step!(a, ra)
        MB.mcs_step!(b, rb)
        a.lattice == b.lattice || return k
    end
    0
end

@testset "The Hamiltonian and T_cpm are jointly identifiable only up to scale" begin

    # ---- the premise, tested before the conclusion --------------------------
    #
    # Bit-exactness holds only if the scaling is genuinely LINEAR PER TERM. Were any
    # term a product of two scaled quantities -- a coefficient times a scaled field --
    # the scaling would be quadratic and the trajectory equality below would fail for
    # a reason this file would otherwise attribute to the acceptance rule. So check
    # f(2θ) == 2·f(θ) exactly, on the terms themselves, first.
    @testset "each term is exactly linear in its own coefficient" begin
        s1, s2 = si_state(SR, 1.0, 1.0), si_state(SR, 2.0, 2.0)
        @test s1.lattice == s2.lattice          # same configuration, different weights
        rng = MersenneTwister(11)
        nadh = nvol = nrad = 0
        # Linearity holds algebraically; the test's value is that it would catch a
        # FUTURE term that is not linear, and it only catches one it exercises. So
        # record what each family was actually sampled over, and require breadth.
        adh_pairs = Set{Tuple{Int,Int}}()
        vol_signs = Set{Int}()
        rad_pairs = Set{Tuple{Int,Int}}()
        sp_of(s, σ) = σ > 0 ? s.cells[Int(σ)].species : 0
        for _ in 1:4000
            sx, sy, sz = rand(rng, 1:SI_P.N), rand(rng, 1:SI_P.N), rand(rng, 1:SI_P.N)
            tx, ty, tz = rand(rng, 1:SI_P.N), rand(rng, 1:SI_P.N), rand(rng, 1:SI_P.N)
            a = SR.compute_delta_H_terms(s1, sx, sy, sz, tx, ty, tz)
            b = SR.compute_delta_H_terms(s2, sx, sy, sz, tx, ty, tz)
            σs, σt = s1.lattice[sx, sy, sz], s1.lattice[tx, ty, tz]
            if a.adh != 0
                nadh += 1; @test b.adh == 2 * a.adh
                push!(adh_pairs, (sp_of(s1, σs), sp_of(s1, σt)))
            end
            if a.vol != 0
                nvol += 1; @test b.vol == 2 * a.vol
                push!(vol_signs, Int(sign(a.vol)))
            end
            if a.rad != 0
                nrad += 1; @test b.rad == 2 * a.rad
                push!(rad_pairs, (sp_of(s1, σs), Int(sign(a.rad))))
            end
            # Exactness of 2x needs the operands in normal range; a subnormal or an
            # overflow would make the equalities above meaningless rather than false.
            @test all(isfinite, (a.adh, a.vol, a.rad, b.adh, b.vol, b.rad))
            @test !any(x -> x != 0 && issubnormal(x), (a.adh, a.vol, a.rad))
        end
        # None of the three assertions above runs if its term is always zero, and a
        # loop that asserts nothing passes. Require each to have actually fired.
        @test nadh > 100 && nvol > 100 && nrad > 100
        # And require breadth, not just volume: all 1436 radiation samples arriving
        # from one species pair would exercise one table entry and read as coverage.
        @test length(adh_pairs) >= 20                 # many distinct species contacts
        @test vol_signs == Set([-1, 1])               # cells both above and below target
        @test length(unique(first.(rad_pairs))) == 8  # medium plus all seven species
    end

    # ---- the trajectory ------------------------------------------------------
    base = si_trace(1.0, 1.0)
    @testset "scaling H and T together reproduces the trajectory bit-exactly" begin
        # A test that runs two simulations and finds them identical proves nothing if
        # neither moved. Require the run to have rearranged the lattice first.
        @test count(base[1] .!= si_trace(1.0, 1.0; n = 0)[1]) > 500
        @test all(iszero, si_state(SR, 1.0, 1.0).melanin)   # ΔH_mel is identically zero here
        @test si_same(base, si_trace(2.0, 2.0))
    end

    # ---- control 1: the asymmetry -------------------------------------------
    #
    # Scaling H without T changes every acceptance probability. If this passed, the
    # test above would be reporting that two identical simulations are identical.
    @testset "control: scaling H without T must diverge" begin
        ctl = si_trace(2.0, 1.0)
        @test !si_same(base, ctl)
        # Not just "they differ" -- the asymmetry changes every acceptance probability,
        # so it must bite in the FIRST sweep. A late divergence would mean something is
        # very nearly cancelling, which would be a different finding.
        @test si_first_divergence(SR, SR, 2.0, 1.0, 1.0, 1.0) == 1
        @test si_first_divergence(SR, SR, 2.0, 2.0, 1.0, 1.0) == 0   # symmetric: never
    end

    # ---- control 2: the unscaled coefficient --------------------------------
    #
    # ΔH_mel is linear in a hard-coded 0.5 literal -- `0.5 * M_local`, currently at
    # biofilms_potts.jl:593 and :596, though grep for the expression rather than
    # trusting those numbers (SCALE-03) -- and not in any parameter, so the symmetric
    # scaling above cannot scale it. With a
    # non-zero melanin field the degeneracy must therefore BREAK. That is a
    # measurement of PP-T2-29 rather than a remark about it, and it makes this test a
    # standing detector: any future coefficient in H that is not reachable through the
    # parameter interface breaks this symmetry and fails here.
    @testset "control: a coefficient outside the parameter interface breaks the symmetry" begin
        mel_base = si_trace(1.0, 1.0; melanin = 1.0)
        mel_scaled = si_trace(2.0, 2.0; melanin = 1.0)
        @test !si_same(mel_base, mel_scaled)
        @test si_first_divergence(SR, SR, 2.0, 2.0, 1.0, 1.0; melanin = 1.0) > 0

        # THE POSITIVE CONTROL, and without it this testset shows only that melanin
        # being present breaks the symmetry -- not that the UNSCALED LITERAL does.
        # Scale the 0.5 as well, in a second copy of the source patched in memory, and
        # the symmetry must return. The substitution is asserted to have applied, and
        # at both sites: a replacement that silently matched nothing would report "no
        # difference" and read exactly like the result being sought.
        @test SI_MEL_SUBS[] == 2
        @test si_first_divergence(SI_SR_MEL2, SR, 2.0, 2.0, 1.0, 1.0; melanin = 1.0) == 0
    end
end
