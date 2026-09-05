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

    @testset "provenance is declared and sticky, not inferred from X_total" begin
        # Codex on 2c93d1e.  Inferring "gated" from X_total != 1.0 fails twice:
        # a coupled state can hit mean(X_tot) == 1.0 by coincidence, and c/s are
        # path-dependent so they stay gated after the value moves back.
        rp1 = SR.RadiolysisParams(Nr = 20, Ddot_R = 1.0, c_ext = 1.0)

        # (a) THE COLLISION: X_total == 1.0, but the basis came from occupancy.
        rd = SR.init_radiolysis(rp1; R = 10.0)
        rd.basis_from_occupancy = true
        @test rd.params.X_total == 1.0        # indistinguishable by value
        @test_throws ErrorException SR.step_radiolysis!(rd, 0.5)

        # (b) STICKINESS: gated, stepped, then the value returns to the default.
        rd2 = SR.init_radiolysis(
            SR.RadiolysisParams(Nr = 20, X_total = 0.065, Ddot_R = 1.0,
                                c_ext = 1.0, basis_gate_ack = true); R = 10.0)
        rd2.basis_from_occupancy = true
        SR.step_radiolysis!(rd2, 0.5)                    # acked, so it runs
        rd2.params = SR.RadiolysisParams(Nr = 20, X_total = 1.0, Ddot_R = 1.0,
                                         c_ext = 1.0)    # ack dropped, value reset
        @test rd2.params.X_total == 1.0
        @test rd2.basis_from_occupancy                   # provenance survives
        @test_throws ErrorException SR.step_radiolysis!(rd2, 0.5)

        # (c) and a genuinely standalone state is still unmarked and allowed.
        rd3 = SR.init_radiolysis(rp1; R = 10.0)
        @test rd3.basis_from_occupancy === false
        SR.step_radiolysis!(rd3, 0.5)
        @test all(isfinite, rd3.c)
    end

    @testset "the coupled path actually sets provenance" begin
        # (a) and (b) above set the flag by hand; this proves the real coupled
        # loop sets it, so the guard is not protecting a field nobody writes.
        params = SR.CPMParams(N = 20, n_cells_per_species = 2,
                              snapshot_interval = 100)
        rp = SR.RadiolysisParams(Nr = 20, Ddot_R = 1.0, c_ext = 1.0,
                                 basis_gate_ack = true)
        sim = SR.init_coupled_simulation(params, rp; seed = 3)
        @test sim.rd.basis_from_occupancy === false   # before any step
        SR.advance_window!(sim, 2)
        @test sim.rd.basis_from_occupancy === true    # after the installer runs
    end

    @testset "the ack census: exactly these sites, no more" begin
        # Codex on 2c93d1e: the first version grepped for the literal
        # "basis_gate_ack = true", so `basis_gate_ack=true`, extra spacing, a
        # line break, or a field assignment all opened the gate while the census
        # stayed green -- a check asserting a bound it could not enforce.
        #
        # This parses instead.  _opens_gate walks the Julia AST for assignments
        # to basis_gate_ack in any syntactic form.
        #
        # SCOPE, stated rather than implied: this is a STATIC check.  A value
        # computed at runtime and splatted in (`RadiolysisParams(; kw...)`) is
        # not statically detectable and this cannot catch it.  What it does
        # establish is that no site opens the gate by writing the field
        # directly, in any spelling.  The known-bad fixtures below prove it
        # catches each form rather than passing because it matches nothing.
        function _opens_gate(ex)
            hit = false
            walk(e) = begin
                if e isa Expr
                    if e.head === :kw || e.head === :(=)
                        lhs, rhs = e.args[1], e.args[2]
                        names = lhs === :basis_gate_ack ||
                                (lhs isa Expr && lhs.head === :. &&
                                 lhs.args[2] == QuoteNode(:basis_gate_ack))
                        # PROPAGATION is not OPENING.  The reconstruction sites
                        # write `basis_gate_ack = rp.basis_gate_ack`, forwarding
                        # a flag someone else already justified; and `= false`
                        # is the closed default. Anything else -- a literal
                        # true, or a computed value whose provenance this cannot
                        # see -- counts as opening the gate and must be listed.
                        forwards = rhs isa Expr && rhs.head === :. &&
                                   rhs.args[2] == QuoteNode(:basis_gate_ack)
                        names && !forwards && rhs !== false && (hit = true)
                    end
                    foreach(walk, e.args)
                end
            end
            walk(ex)
            return hit
        end
        opens(src) = _opens_gate(Meta.parseall(src))

        # The control: each of these bypassed the old literal grep.
        @test opens("RadiolysisParams(basis_gate_ack=true)")
        @test opens("RadiolysisParams(basis_gate_ack  =  true)")
        @test opens("RadiolysisParams(Nr = 20,\n  basis_gate_ack =\n  true)")
        @test opens("rd.params.basis_gate_ack = true")
        @test opens("x = (basis_gate_ack = true,)")
        # and it must not fire on the guard's own comparison or on prose,
        # or the census would name every file that mentions the symbol.
        @test !opens("ack === :determinism_only")
        @test !opens("f(basis_gate_ack)")
        @test !opens("# basis_gate_ack = true, in a comment")
        @test !opens("s = \"basis_gate_ack = true\"")
        # propagation forwards a flag already justified elsewhere; the closed
        # default opens nothing. Neither should name a site.
        @test !opens("RadiolysisParams(basis_gate_ack = rp.basis_gate_ack)")
        @test !opens("basis_gate_ack::Bool = false")
        # but a value this cannot trace IS an opening, and must be listed
        @test opens("RadiolysisParams(basis_gate_ack = should_ack())")

        expected = Set([
            "validate_serial.jl",              # CPM trajectory determinism
            "tests/genealogy_tests.jl",        # legacy vs windowed API equivalence
            "tests/checkpoint_io_tests.jl",    # snapshot/restart round trip
            "export_checkpoint.jl",            # interchange export; labels the file
            "tests/radiodialysis_basis_gate.jl",  # this file, testing the gate
            "tests/jacc_parity_tests.jl",      # CPM acceptance counts; the
                                              # 10x-uptake testset in that file
                                              # holds the claim to account
            "jacc_acceptance_figure.jl",      # fig5; plots the same acceptance
                                              # quantities, same claim
            "regenerate_fig3.jl",             # fig3 is m(t) and P_eff/P0, which
                                              # the gated basis does not reach;
                                              # that script proves it before it
                                              # copies anything
        ])
        repo = dirname(@__DIR__)
        found = Set{String}()
        for (root, _, files) in walkdir(repo)
            occursin("/.git", root) && continue
            for fname in files
                endswith(fname, ".jl") || continue
                path = joinpath(root, fname)
                src = read(path, String)
                parsed = try
                    Meta.parseall(src)
                catch
                    continue
                end
                _opens_gate(parsed) && push!(found, relpath(path, repo))
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
