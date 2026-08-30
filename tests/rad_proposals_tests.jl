# ---------- The producer for §6.2's per-proposal ΔH_rad statistics, and the
# ---------- proof that instrumenting them does not move the trajectory.
#
# §6.2 states 26.7% of evaluated proposals carry ΔH_rad < -5e-5 over 3968838
# proposals across seeds 42, 43 and 44 at 400 MCS. Version 1.2 stated 29.8% over
# an unnamed run and nothing in the repository could produce either figure. This
# pins the restated one to the code that makes it.

include(joinpath(REPO, "rad_proposals.jl"))

@testset "§6.2's proposal statistics reproduce from the shipped hook" begin
    SR_rp = load_serial()

    @testset "the hook does not move the trajectory" begin
        # SAME SHAPE AS `the driver tally is inert`. `on_proposal` is called
        # before the acceptance draw, so a harness that consulted the generator
        # would reorder every subsequent move; the shipped one does not, and
        # this is what says so.
        function run(hooked::Bool)
            p = SR_rp.CPMParams(N = 16, n_cells_per_species = 3)
            st = SR_rp.init_state(p; seed = 5)
            rng = Random.Xoshiro(5)
            seen = Ref(0)
            for _ in 1:6
                SR_rp.mcs_step!(st, rng;
                    on_proposal = hooked ? ((t, dh) -> (seen[] += 1)) : nothing)
                SR_rp.update_melanin!(st)
            end
            return copy(st.lattice), copy(st.melanin), seen[]
        end
        bare_lat, bare_mel, bare_seen = run(false)
        hook_lat, hook_mel, hook_seen = run(true)
        @test bare_lat == hook_lat
        @test bare_mel == hook_mel

        # THE HOOK MUST HAVE FIRED. Equality between a bare run and a run whose
        # hook was never called is vacuous, and that is how a broken hook looks.
        @test bare_seen == 0
        @test hook_seen > 0

        # CONTROL 1, AND WHAT IT IS FOR: a different seed must give a different
        # lattice. This establishes that the COMPARISON varies with input at all,
        # so equality at seed 5 is not vacuously true. It says nothing about
        # whether the hook channel could have disturbed anything.
        other = SR_rp.init_state(SR_rp.CPMParams(N = 16, n_cells_per_species = 3); seed = 6)
        rng6 = Random.Xoshiro(6)
        for _ in 1:6
            SR_rp.mcs_step!(other, rng6)
            SR_rp.update_melanin!(other)
        end
        @test other.lattice != bare_lat

        # CONTROL 2, WHICH IS THE ONE THAT MAKES THE NULL READABLE. Firing and
        # firing WITH THE POWER TO DISTURB are different properties, and only the
        # second makes the equality above evidence about `on_proposal`. The JACC
        # port's guarantee was structural -- a counter-based RNG keyed on (seed,
        # step, site) has no stream position to advance -- and this file does not
        # inherit it: `mcs_step!` draws from a sequential generator and the hook
        # sits immediately before the acceptance draw. So the claim is empirical
        # and needs its own sensitivity evidence.
        #
        # BOTH READINGS OF A NULL, DECIDED BEFORE IT IS SEEN. If the drawing hook
        # does NOT diverge, that is ambiguous between "the comparison is blind"
        # and "`on_proposal` has no reach into the stream", which would be the
        # stronger result. Object identity separates them, so it is asserted
        # first: `mcs_step!` consumes from the generator it is passed and
        # re-seeds nothing, therefore divergence is the expected outcome and a
        # null would mean the harness, not the hook.
        rng5 = Random.Xoshiro(5)
        drew = SR_rp.init_state(SR_rp.CPMParams(N = 16, n_cells_per_species = 3); seed = 5)
        captured = Ref{Any}(nothing)
        greedy = (t, dh) -> (captured[] = rng5; rand(rng5); nothing)
        for _ in 1:6
            SR_rp.mcs_step!(drew, rng5; on_proposal = greedy)
            SR_rp.update_melanin!(drew)
        end
        @test captured[] === rng5          # same object the stepper draws from
        @test drew.lattice != bare_lat     # the hook channel CAN perturb
    end

    @testset "the published numbers" begin
        tn = 0; tb = 0; lo = Inf; hi = -Inf
        for seed in (42, 43, 44)
            r = Base.invokelatest(rad_proposals, SR_rp, seed, 400)
            ref = PUBLISHED[seed]
            @test r.n == ref.n
            @test r.below == ref.below
            @test r.min ≈ ref.min
            tn += r.n; tb += r.below; lo = min(lo, r.min); hi = max(hi, r.max)
        end
        @test tn == POOLED_N
        @test round(tb / tn; digits = 3) == POOLED_FRAC

        # The acceptance-favouring extreme at N=40 is the single-role value, and
        # the pairwise extremum is reached only in the DISFAVOURING direction.
        # That asymmetry is the sentence §6.2 makes, so it is pinned here: if
        # `lo` ever reached -0.07505 the manuscript would be wrong.
        @test lo ≈ -0.075
        @test hi ≈ 0.07505
        @test lo > -0.07505
    end

    @testset "the manuscript states what this measures" begin
        tex = read(joinpath(REPO, "preprint",
                   "modeling_radioresistance_and_radiotropic_fitness.tex"), String)
        @test occursin("26.7", tex)
        @test occursin(raw"3\,968\,838", tex)

        # THE DEFECT WAS THE UNNAMED RUN, NOT THE DIGITS. A correction record has
        # to be able to name the number it withdraws -- 29.8% still appears once,
        # in the version 1.2 note, and forbidding the string outright would forbid
        # the paper from recording its own correction. What must not survive is
        # the phrasing that reported a measurement over a run nothing identified.
        @test !occursin("Measured over the run", tex)
        @test occursin(raw"Measured over the three runs of Table~\ref{tab:decided_moves}", tex)
        @test count(_ -> true, eachmatch(r"29\.8", tex)) == 1
    end
end
