# ---------- The prose gate: a bound stated in the paper must be the bound the
# ---------- code has.
#
# Section 6.2 asserted that the largest acceptance-favouring contribution of
# ΔH_rad was 5e-5, "so ... reversing the move would require the drawn variate to
# fall within about 1e-5 of the acceptance threshold", and concluded from that
# that observing zero radiation-decided moves was analytically expected rather
# than an under-sampling artifact. The bound was wrong by 1.5e3 and no check in
# the repository could have said so: contract_csv.jl guards the trajectory,
# delta_h_decomposition.jl guards that the four terms sum to the scalar, and
# test_claims_ledger.py guards that retracted PHRASES do not survive. None of
# them reads a NUMBER out of the prose and compares it to the model.
#
# WHAT MAKES THE BOUND CHECKABLE. ΔH_rad is +β_ion[source]·I when the source
# gains a site and -β_ion[target]·I when the target loses one, so a move is
# favoured whenever the total is negative and BOTH SIGNS OF β_ion CAN DO THAT --
# a positively signed species vacating a site favours acceptance exactly as much
# as a negatively signed one occupying it. The bound is therefore
# max|β_ion| · I₀ over species, which is a property of the shipped coefficients
# and not of any run. The previous bound took the minimum magnitude among the
# negatively signed species, which is the same arithmetic with one role missing.
#
# SCOPE, STATED RATHER THAN IMPLIED. This gates ONE number: the ΔH_rad bound in
# §6.2. It is not a general prose-versus-code checker and does not pretend to
# be. What it establishes is that this particular class of defect -- a bound
# asserted in prose that the coefficients contradict -- now has somewhere to
# fail, for the one bound that has already failed once.

const TEX = joinpath(REPO, "preprint",
                     "modeling_radioresistance_and_radiotropic_fitness.tex")

"Parse a LaTeX scientific literal, e.g. `7.5\\times10^{-2}` -> 0.075."
function _latex_number(s::AbstractString)
    m = match(r"([0-9]*\.?[0-9]+)\s*\\times\s*10\^\{(-?[0-9]+)\}", s)
    m === nothing && return nothing
    return parse(Float64, m[1]) * 10.0^parse(Int, m[2])
end

"The bound §6.2 states for the acceptance-favouring reach of ΔH_rad."
function stated_rad_bound(tex::AbstractString)
    m = match(r"\\max_s\s*\|\\beta_\{s,\\mathrm\{ion\}\}\|\s*=\s*([^,\$]+)", tex)
    m === nothing && return nothing
    return _latex_number(m[1])
end

@testset "the ΔH_rad bound in prose is the bound in the coefficients" begin
    @test isfile(TEX)
    tex = read(TEX, String)
    p = SR.CPMParams()

    # BOTH ROLES, WHICH IS THE WHOLE CORRECTION. Not `minimum(abs, ...)` over the
    # negatively signed entries, which is what the withdrawn sentence used.
    true_bound = maximum(abs, p.β_ion) * p.I0
    stated = stated_rad_bound(tex)

    @test stated !== nothing
    @test isapprox(stated, true_bound; rtol = 1e-6)

    # ...and the same sentence must not still be asserting the withdrawn one.
    @test !occursin("largest radiation contribution that could favour an accepted move is of that order", tex)

    @testset "the control: the version 1.1 sentence must fail this" begin
        # Committed inline rather than recovered with `git show`: a squash-merge
        # makes a pinned sha unreachable and the control degrades into a skip.
        withdrawn = raw"""
        $\beta_{s,\mathrm{ion}}$ is negative for two of the seven species at $-5\times10^{-5}$, so
        the largest radiation contribution that could favour an accepted move is of that order,
        and reversing the move would require the drawn variate to fall within about $10^{-5}$ of
        the acceptance threshold."""
        # It states no max| | at all, so the extractor returns nothing -- which
        # is a failure of this gate, not a pass. Assert that explicitly, or
        # "no bound stated" and "the right bound stated" would look identical.
        @test stated_rad_bound(withdrawn) === nothing

        # And a sentence that DOES state the withdrawn magnitude must not agree
        # with the coefficients.
        wrong = raw"$\max_s |\beta_{s,\mathrm{ion}}| = 5\times10^{-5}$"
        @test stated_rad_bound(wrong) ≈ 5e-5
        @test !isapprox(stated_rad_bound(wrong), true_bound; rtol = 1e-6)
    end

    @testset "the parser reads what LaTeX actually writes" begin
        # isapprox, not ==: 7.5 * 10.0^-2 is 0.07500000000000001, and pinning
        # the parser to exact binary equality would fail on arithmetic rather
        # than on parsing.
        @test _latex_number(raw"7.5\times10^{-2}") ≈ 7.5e-2
        @test _latex_number(raw"5\times10^{-5}")   ≈ 5e-5
        @test _latex_number(raw"1.5\times10^{3}")  ≈ 1.5e3
        @test _latex_number("no number here")      === nothing
    end
end
