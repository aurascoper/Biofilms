# The RADIODIALYSIS: BLOCKED basis gate, its one exemption, and the R/Julia
# asymmetry it pins.
#
# WHY THIS IS AN ASYMMETRY TEST AND NOT AN AGREEMENT TEST.  The two
# implementations deliberately DISAGREE at non-unit X_total, and that has to be
# asserted rather than left to be discovered:
#
#   biofilms_radiodialysis.R  derives X_red = f_red_active * X_total
#                             (uptake_rate_of), so U scales with the basis and
#                             lambda ~ 1/sqrt(U).  Asserted there by check 15 of
#                             analysis/verify_biofilm_depth_profile.R.
#   this path                 REFUSES, because the X_red reaching it is one
#                             species' occupied sites over ALL interior sites --
#                             neither a biomass fraction nor a reducer fraction
#                             (README.md:344), which is the defect the gate names.
#
# Making this path conformant without repairing that quantity would turn a
# visible defect into an invisible one.  So the refusal is asserted here, and
# this test fails the day it is removed -- forcing whoever removes it to confront
# the quantity underneath rather than only the arithmetic on top.
#
# Before this gate existed the block was an ACCIDENT: the additive uptake form
# happens to misbehave at non-unit X_total, and that side effect was standing in
# for a guard.  A tripwire that any tidy-up can remove is not a gate.

@testset "radiodialysis basis gate" begin

    @testset "the default basis still integrates" begin
        # Without this the refusal half is vacuous: a gate that refused
        # unconditionally would satisfy every assertion below it.
        rp = SR.RadiolysisParams(Nr = 20, Ddot_R = 1.0, c_ext = 1.0)
        @test rp.X_total == 1.0
        @test rp.basis_gate_ack === false
        rd = SR.init_radiolysis(rp; R = 10.0)
        SR.step_radiolysis!(rd, 0.5)
        @test all(isfinite, rd.c)
    end

    @testset "a coupled basis is refused, by name" begin
        for xt in (0.065, 0.1, 0.5, 0.9, 1.5, 2.0)
            rp = SR.RadiolysisParams(Nr = 20, X_total = xt,
                                     Ddot_R = 1.0, c_ext = 1.0)
            rd = SR.init_radiolysis(rp; R = 10.0)
            err = try
                SR.step_radiolysis!(rd, 0.5); nothing
            catch e
                sprint(showerror, e)
            end
            @test err !== nothing
            @test occursin("RADIODIALYSIS: BLOCKED", err)
            @test occursin("neither a biomass fraction nor a reducer fraction", err)
        end
    end

    @testset "the refusal precedes any stepping" begin
        # If the arithmetic were still doing the blocking, state would move
        # before anything threw.
        rp = SR.RadiolysisParams(Nr = 20, X_total = 0.065,
                                 Ddot_R = 1.0, c_ext = 1.0)
        rd = SR.init_radiolysis(rp; R = 10.0)
        c0, s0, t0 = copy(rd.c), copy(rd.s), rd.t
        @test_throws ErrorException SR.step_radiolysis!(rd, 0.5)
        @test rd.c == c0
        @test rd.s == s0
        @test rd.t == t0
    end

    @testset "the determinism exemption is narrow" begin
        rp = SR.RadiolysisParams(Nr = 20, X_total = 0.065, Ddot_R = 1.0,
                                 c_ext = 1.0, basis_gate_ack = true)
        rd = SR.init_radiolysis(rp; R = 10.0)
        SR.step_radiolysis!(rd, 0.5)          # allowed
        @test all(isfinite, rd.c)
        # and ONLY that symbol opens it
        # and the default still refuses
        rpx = SR.RadiolysisParams(Nr = 20, X_total = 0.065, Ddot_R = 1.0,
                                  c_ext = 1.0)
        @test rpx.basis_gate_ack === false
        rdx = SR.init_radiolysis(rpx; R = 10.0)
        @test_throws ErrorException SR.step_radiolysis!(rdx, 0.5)
    end

    @testset "the ack census: exactly these sites, no more" begin
        # An exemption that can be added silently is not an exemption, it is a
        # default.  Every site opening the gate is enumerated here, so widening
        # it is a deliberate edit to this list rather than a line someone adds.
        expected = Set([
            "validate_serial.jl",              # CPM trajectory determinism
            "tests/genealogy_tests.jl",        # legacy vs windowed API equivalence
            "tests/checkpoint_io_tests.jl",    # snapshot/restart round trip
            "tests/radiodialysis_basis_gate.jl",  # this file, testing the gate
        ])
        repo = dirname(@__DIR__)
        found = Set{String}()
        for (root, _, files) in walkdir(repo)
            occursin("/.git", root) && continue
            for fname in files
                endswith(fname, ".jl") || continue
                path = joinpath(root, fname)
                # The ASSIGNMENT form only.  The guard itself compares with
                # `ack === :determinism_only`, and the docstrings name the
                # symbol in prose; neither opens the gate for anyone.
                if occursin("basis_gate_ack = true", read(path, String))
                    push!(found, relpath(path, repo))
                end
            end
        end
        @test found == expected
    end

    @testset "the exemption's claim is true: nothing recorded depends on the basis" begin
        # validate_serial.jl acknowledges a blocked basis on the grounds that it
        # records no radiodialysis quantity.  That is a CLAIM, and this is what
        # holds it to account: the gated basis enters only through `uptake`, so
        # changing the uptake constants by 10x must leave every recorded column
        # byte-identical.  The day someone records a c- or s-derived quantity,
        # this fails and the exemption has to be re-argued.
        function recorded(k_ads, k_red)
            params = SR.CPMParams(N = 40, n_cells_per_species = 6,
                                  snapshot_interval = 20)
            rp = SR.RadiolysisParams(Nr = 40, Ddot_R = 1.0, c_ext = 1.0,
                                     k_ads = k_ads, k_red = k_red,
                                     basis_gate_ack = true)
            state, rd, _, _ = SR.run_simulation_coupled(params, rp, 40; seed = 42)
            snap = SR.take_snapshot(state, 40)
            rows = [(sd.species, sd.total_volume, sd.n_cells, sd.mean_melanin)
                    for sd in snap.species_data]
            return (rows, length(state.cells), rd.m)
        end
        base = recorded(0.05, 0.02)
        perturbed = recorded(0.5, 0.2)        # 10x the uptake constants
        @test base == perturbed
    end
end
