# Console report honesty: no fabricated physical units for uncalibrated
# quantities. PR #12 finding: the banner honestly admitted Ddot_R has no
# seconds/MCS conversion, but the "Final P_eff" and "Cumulative dose at
# membrane" lines a few lines later still printed cm/s / Gy as if they were
# real measurements. `print_membrane_report` is split out of `main_coupled()`
# specifically so it's testable without a full 100-MCS simulation run (and
# without `export_figures`, which `main_coupled()` calls unconditionally and
# which lives below the CairoMakie split marker, out of reach of the test
# sandbox module `SR`).
#
# `rd.t` is set nonzero so the printed numbers are non-trivial, not just the
# zero-time initial condition.

let
    rp = SR.RadiolysisParams(Nr = 10, Ddot_R = 1.0, c_ext = 1.0)
    rd = SR.init_radiolysis(rp)
    rd.t = 50.0

    old_stdout = stdout
    pipe_rd, pipe_wr = redirect_stdout()
    try
        SR.print_membrane_report(rp, rd, 100)
    finally
        redirect_stdout(old_stdout)
        close(pipe_wr)
    end
    report = String(read(pipe_rd))
    close(pipe_rd)

    # THE DETECTORS FIRST. An absence assertion run only against corrected
    # output is a check that cannot fail: weaken either pattern -- a typo, a
    # stray escape, `cm/s` drifting to `cm s^-1` -- and it stops matching
    # anything at all, while the test stays green and reports the fabricated
    # units are gone. So give each one the known-bad string it was written to
    # catch, taken from the report as it actually printed before the fix, and
    # require it to bite before its absence from production output means
    # anything.
    P_EFF_UNITS = r"Final P_eff = [\d.]+ cm/s"
    DOSE_UNITS  = r"Cumulative dose at membrane: [\d.]+ Gy"

    known_bad_p_eff = "  Final P_eff = 0.0123 cm/s\n"
    known_bad_dose  = "  Cumulative dose at membrane: 4.56 Gy\n"
    @test occursin(P_EFF_UNITS, known_bad_p_eff)
    @test occursin(DOSE_UNITS, known_bad_dose)
    # And each must be specific to its own line, not a pattern loose enough to
    # fire on any report text -- which would make the absence assertions below
    # pass for the wrong reason too.
    @test !occursin(P_EFF_UNITS, known_bad_dose)
    @test !occursin(DOSE_UNITS, known_bad_p_eff)

    # Only now: the two specific fabricated lines this finding named are gone
    # from what `print_membrane_report` actually prints.
    @test !occursin(P_EFF_UNITS, report)
    @test !occursin(DOSE_UNITS, report)

    # The honest replacements must be present.
    @test occursin("not Gy", report)
    @test occursin("not computed in this work", report)
end

# ---------------------------------------------------------------- every producer
#
# THE FIRST FIX COVERED ONE HELPER AND THE SUITE ONLY EXERCISED THAT HELPER, so
# the summary honestly declined to state an absolute permeability while two live
# `[RD]` lines in `main_coupled` and the Shiny status in
# `biofilms_radiodialysis.R` kept printing one -- the R line with `cm/s`
# attached, which is the fabricated unit the whole finding was about. A reader
# running either program still received the withdrawn value through an ordinary
# report path.
#
# SCAN EVERY LINE, NOT ONLY LINES CONTAINING `printf`. The first version of this
# scan filtered to lines mentioning printf/sprintf, and in R the format string
# sits on its own line below the `sprintf(` call -- so the offending line was
# never examined and the negative control passed against the reverted bug. A
# filter that narrows the search to where the author expected the defect is not
# a search.

let
    sources = [joinpath(REPO, f) for f in
               ("biofilms_potts.jl", "biofilms_radiodialysis.R")]

    # `P_eff=` / `P_eff(t)=` followed by a format specifier: an absolute value
    # about to be printed. The ratio forms write `P_eff/P0=` and do not match.
    ABSOLUTE = r"P_eff(\([^)]*\))?\s*=\s*%"

    # The detector must find a planted hit before its absence means anything.
    @test occursin(ABSOLUTE, "  \"P_eff=%.5f cm/s\",")
    @test occursin(ABSOLUTE, "  \"P_eff(t_end)=%.4f cm/s\",")
    @test !occursin(ABSOLUTE, "  \"P_eff/P0=%.5f\",")
    @test !occursin(ABSOLUTE, "  \"P_eff(t_end)/P0=%.4f (dimensionless)\",")

    offenders = Tuple{String,Int,String}[]
    for path in sources
        for (i, line) in enumerate(eachline(path))
            startswith(strip(line), "#") && continue
            occursin(ABSOLUTE, line) || continue
            # Allowed ONLY where the line says the absolute value is not
            # computed -- the honest replacement the summary already carries.
            occursin("not computed", line) && continue
            push!(offenders, (basename(path), i, strip(line)))
        end
    end
    if !isempty(offenders)
        for (f, i, l) in offenders
            println("  absolute P_eff still reported: $f:$i  $(first(l, 70))")
        end
    end
    @test isempty(offenders)

    # And no live report may attach cm/s to a permeability number.
    unit_hits = Tuple{String,Int}[]
    for path in sources
        for (i, line) in enumerate(eachline(path))
            startswith(strip(line), "#") && continue
            (occursin("P_eff", line) && occursin("cm/s", line) &&
             !occursin("not computed", line)) || continue
            push!(unit_hits, (basename(path), i))
        end
    end
    @test isempty(unit_hits)
end
