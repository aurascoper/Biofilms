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

    # The two specific fabricated lines this finding named must be gone.
    @test !occursin(r"Final P_eff = [\d.]+ cm/s", report)
    @test !occursin(r"Cumulative dose at membrane: [\d.]+ Gy", report)

    # The honest replacements must be present.
    @test occursin("not Gy", report)
    @test occursin("not computed in this work", report)
end
