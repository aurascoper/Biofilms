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
# as a negatively signed one occupying it.
#
# THE BOUND IS OVER PAIRINGS, NOT OVER SPECIES, which is the second thing this
# gate got wrong. A copy between two occupied parcels contributes
# (β_source - β_target)·I, so the acceptance-favouring reach is the extremum over
# the (source, target) pairings the lattice permits, not the largest single
# coefficient. max|β_ion|·I₀ = 7.5e-2 takes one role at a time and is therefore
# not a bound at all: a CS source (-5e-5) copying into an SO target (7.5e-2)
# reaches -7.505e-2, and DOES SO IN A SHIPPED CONFIGURATION -- four times in
# 1298668 evaluated proposals at N=20, seed 42, 400 MCS. Raised by Codex on pull
# request #23; the withdrawn version of this file computed max|β_ion| and would
# have rejected the correct value.
#
# max(0, ·) ON EACH ROLE rather than (max β - min β): an absent parcel
# contributes nothing, so if every coefficient shared one sign the bound would
# come from one role alone. Writing the difference directly would bake in a
# property of today's seven coefficients that no test states.
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
    # THE WHOLE EXPRESSION, NOT A LANDMARK AND A NUMBER. The first version of
    # this required only that \max_{s,t} be followed somewhere by \beta_{t,...}
    # and the right literal, so `\max_{s,t}\beta_t = 7.505e-2` and
    # `\max_{s,t}(\beta_s + \beta_t)I_\gamma = 7.505e-2` both passed every
    # assertion below. A semantic regression could keep the corrected number and
    # invalidate the bound it names, which is the failure this file exists to
    # catch. Raised by Codex on pull request #23.
    #
    # So the shape is required in full: both subscripts in their roles, the
    # subtraction between them in that order, and I_γ. The single-species form
    # stays unmatched on purpose -- a sentence bounding max_s |β| is one role at
    # a time and reading a number out of it would launder that.
    m = match(r"""\\max_\{s,t\}\s*\(\s*\\beta_\{t,\\mathrm\{ion\}\}\s*-\s*\\beta_\{s,\\mathrm\{ion\}\}\s*\)\s*(?:\\,)?\s*I_\\gamma\s*=\s*([^,\$]+)"""x, tex)
    m === nothing && return nothing
    return _latex_number(m[1])
end

@testset "the ΔH_rad bound in prose is the bound in the coefficients" begin
    @test isfile(TEX)
    tex = read(TEX, String)
    p = SR.CPMParams()

    # BOTH ROLES AND BOTH ENDS OF THE PAIRING. Not `minimum(abs, ...)` over the
    # negatively signed entries (version 1.1), and not `maximum(abs, ...)` over
    # species (version 1.2): see the header.
    βmax, βmin = maximum(p.β_ion), minimum(p.β_ion)
    true_bound = (max(0.0, βmax) + max(0.0, -βmin)) * p.I0
    stated = stated_rad_bound(tex)

    @test stated !== nothing
    @test isapprox(stated, true_bound; rtol = 1e-6)

    # THE BOUND IS ATTAINED, NOT MERELY ASSERTED. Both ends come from a distinct
    # species, so two distinct parcels can realize it; if argmax and argmin
    # coincided the extremum would need one parcel in both roles at once, which
    # the copy rule forbids, and the bound would be loose without saying so.
    @test argmax(p.β_ion) != argmin(p.β_ion)

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

        # A sentence in the withdrawn SHAPE states no bound on ΔH_rad at all --
        # max_s |β| is one role at a time -- so the extractor must refuse it
        # whatever number it carries. Both the 1.1 magnitude and the 1.2 one.
        @test stated_rad_bound(raw"$\max_s |\beta_{s,\mathrm{ion}}| = 5\times10^{-5}$") === nothing
        @test stated_rad_bound(raw"$\max_s |\beta_{s,\mathrm{ion}}| = 7.5\times10^{-2}$") === nothing
    end

    @testset "the control: pairwise-SHAPED formulas that are not the bound" begin
        # SCOPE: the four synthetic sentences listed in this block, and no
        # location in the manuscript. Manuscript coverage is the dH_sites loop.
        # EACH OF THESE CARRIES THE CORRECT LITERAL AND STATES A DIFFERENT
        # QUANTITY. They are what a semantic regression looks like when the
        # number survives the edit, and the previous extractor accepted all of
        # them. Codex named the first two.
        @test stated_rad_bound(raw"$\max_{s,t}\beta_{t,\mathrm{ion}} = 7.505\times10^{-2}$") === nothing
        @test stated_rad_bound(raw"$\max_{s,t}(\beta_{t,\mathrm{ion}} + \beta_{s,\mathrm{ion}})\,I_\gamma = 7.505\times10^{-2}$") === nothing
        # roles swapped: this is the DISfavouring extremum, a different quantity
        @test stated_rad_bound(raw"$\max_{s,t}(\beta_{s,\mathrm{ion}} - \beta_{t,\mathrm{ion}})\,I_\gamma = 7.505\times10^{-2}$") === nothing
        # I_gamma dropped: a bound on the coefficients, not on ΔH_rad
        @test stated_rad_bound(raw"$\max_{s,t}(\beta_{t,\mathrm{ion}} - \beta_{s,\mathrm{ion}}) = 7.505\times10^{-2}$") === nothing

        # AND THE POSITIVE CONTROL, or the four above pass by the extractor
        # being broken rather than by being strict.
        @test stated_rad_bound(raw"$\max_{s,t}(\beta_{t,\mathrm{ion}} - \beta_{s,\mathrm{ion}})\,I_\gamma = 7.505\times10^{-2}$") ≈ 7.505e-2
    end

    @testset "the control: the version 1.2 bound must fail this" begin
        # 7.505e-2 against 7.5e-2 is 0.067% -- larger than the rtol above, so the
        # assertion can tell them apart. Stated, because a control that cannot
        # separate the two values it names proves nothing about either.
        v12 = maximum(abs, p.β_ion) * p.I0
        @test v12 ≈ 7.5e-2
        @test !isapprox(v12, true_bound; rtol = 1e-6)
        @test true_bound > v12

        wrong12 = raw"$\max_{s,t}(\beta_{t,\mathrm{ion}} - \beta_{s,\mathrm{ion}})\,I_\gamma = 7.5\times10^{-2}$"
        @test stated_rad_bound(wrong12) ≈ 7.5e-2
        @test !isapprox(stated_rad_bound(wrong12), true_bound; rtol = 1e-6)
    end

    @testset "the parser reads what LaTeX actually writes" begin
        # isapprox, not ==: 7.5 * 10.0^-2 is 0.07500000000000001, and pinning
        # the parser to exact binary equality would fail on arithmetic rather
        # than on parsing.
        @test _latex_number(raw"7.505\times10^{-2}") ≈ 7.505e-2
        @test _latex_number(raw"7.5\times10^{-2}") ≈ 7.5e-2
        @test _latex_number(raw"5\times10^{-5}")   ≈ 5e-5
        @test _latex_number(raw"1.5\times10^{3}")  ≈ 1.5e3
        @test _latex_number("no number here")      === nothing
    end
end

# ---------------------------------------------- the melanin/radiation RATIO
#
# A DIFFERENT DEFECT SIGNATURE FROM THE BOUND ABOVE. §6.2 said the corrected
# reach is "only fifteen times smaller than the melanin term rather than four
# orders". Neither ratio the shipped values produce is fifteen: 15.41 is the
# acceptance-bias ratio computed against a reach of 5e-2, a magnitude appearing
# nowhere in the paper. The sentence was written mid-correction -- the "four
# orders" half was already fixed while the ratio still used the pre-pairwise
# number -- so this is not a place a correction failed to REACH, it is a place
# two corrections landed at different times and the earlier input survived
# inside the later one's sentence.
#
# You do not find that by asking where a number went. You find it by
# RECOMPUTING the number from values in the same document, which is what the
# bound check above already does, which is why this belongs beside it.
#
# TWO RATIOS, 6% APART, AND EACH LOCATION USES THE ONE ITS REGISTER DEMANDS.
# §6.2, §7.1 and the Conclusion compare Hamiltonian terms, so they carry the ΔH
# ratio. The abstract quotes 15.5%, an acceptance bias, so it carries the bias
# ratio. An unlabelled number in two places is how the next sweep finds one
# claim disagreeing with itself.

"The melanin coupling is hard-coded in compute_delta_H_terms, not tabulated."
function _melanin_coupling(src::AbstractString)
    m = match(r"ΔH_mel\s*-=\s*([0-9.]+)\s*\*\s*M_local", src)
    m === nothing && return nothing
    return parse(Float64, m[1])
end

@testset "the melanin/radiation ratio in prose is the ratio in the constants" begin
    tex = read(TEX, String)
    p = SR.CPMParams()
    src = read(joinpath(REPO, "biofilms_potts.jl"), String)

    coupling = _melanin_coupling(src)
    @test coupling !== nothing
    @test coupling ≈ 0.5

    M_REPORTED = 1.44          # the melanin level Table 3's row is stated at
    dH_mel = coupling * M_REPORTED
    @test dH_mel ≈ 0.720

    βmax, βmin = maximum(p.β_ion), minimum(p.β_ion)
    reach = (max(0.0, βmax) + max(0.0, -βmin)) * p.I0

    ratio_dH   = dH_mel / reach
    ratio_bias = (exp(dH_mel / p.T_cpm) - 1) / (exp(reach / p.T_cpm) - 1)
    @test isapprox(ratio_dH,   9.6;  atol = 0.05)
    @test isapprox(ratio_bias, 10.2; atol = 0.05)

    # ...and the two are NOT interchangeable, which is why each location names
    # its register. If they ever converge this assertion says so.
    @test !isapprox(ratio_dH, ratio_bias; rtol = 1e-3)

    # EVERY LOCATION, NOT THE FIRST MATCH. `match` returns one hit, so the first
    # version of this read §6.2 alone while the docstring claimed each location
    # carries the ratio its register demands. Mutating §7.1's 9.6 to 15.0 left
    # stated_dH at §6.2's value and the whole suite passed; the abstract and the
    # Conclusion were never inspected at all. A test contract that exceeds its
    # assertion is the F_s defect living in a comment. Raised by Codex on #23.
    #
    # Each location gets its own pattern AND the count is pinned, so a location
    # losing its phrasing fails rather than silently dropping out of the sweep.
    dH_sites = [
        ("6.2", r"smaller than the melanin term by a factor of \$([0-9.]+)\$ in"),
        ("7.1", r"radiation term by about an order of magnitude in \$\\Delta H\$ \(\$([0-9.]+)\$\)"),
    ]
    for (where, pat) in dH_sites
        hits = collect(eachmatch(pat, tex))
        @test length(hits) == 1              # the location exists, exactly once
        for h in hits
            @test isapprox(parse(Float64, h[1]), ratio_dH; atol = 0.05)
        end
    end

    # The Conclusion and the abstract state the ratio in words rather than
    # digits, so they are pinned by phrase. Naming that difference is the point:
    # an unlabelled sweep would report them as covered.
    @test occursin("exceeding the direct\nspecies-specific radiation term by about an order of magnitude in \$\\Delta H\$", tex)
    @test occursin("an acceptance bias about ten\ntimes smaller rather than four orders smaller", tex)

    @testset "the control generalises past the value it was built from" begin
        # A synthetic sentence stating FIFTEEN failing would only re-test the
        # training case. The control is a DIFFERENT wrong ratio, plus the
        # correct one passing -- that pins an arithmetic relation instead of one
        # string, and the next half-correction will not land on fifteen either.
        wrong = raw"smaller than the melanin term by a factor of $12.0$ in $\Delta H$"
        parsed = match(r"smaller than the melanin term by a factor of \$([0-9.]+)\$ in", wrong)
        @test parsed !== nothing
        @test !isapprox(parse(Float64, parsed[1]), ratio_dH; atol = 0.05)

        right = raw"smaller than the melanin term by a factor of $9.6$ in $\Delta H$"
        ok = match(r"smaller than the melanin term by a factor of \$([0-9.]+)\$ in", right)
        @test isapprox(parse(Float64, ok[1]), ratio_dH; atol = 0.05)

        # And fifteen, as a special case rather than as the test.
        @test !isapprox(15.41, ratio_bias; atol = 0.05)
    end
end
