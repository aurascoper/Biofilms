# The manuscript makes several checkable claims about its own source code and its own
# bibliography. Each is currently true only because nobody has changed the code or the
# citations since it was written — nothing enforces it staying true. This file makes four
# of those claims fail loudly the moment they stop holding, instead of silently.
#
# SCOPE. "Every .jl and .R file" in the manuscript's own prose is broader than what any
# assertion here checks, and broader than what the claims are actually about: the three R
# codes and the CPM the manuscript discusses are the repo-root simulation sources, not
# tests/, coupling/ or calibration/ infrastructure. The corrected .tex prose (see
# preprint/...tex, section 3.2) now names this same file set explicitly, so the prose and
# the check assert the same thing -- the original defect was a prose claim broader than
# what was actually enforced, and a mismatched scope here would repeat it.
SIM_FILE_NAMES = [
    "biofilms_potts.jl", "biofilms_potts_jacc.jl", "biofilms.R", "biofilms_3d.R",
    "biofilms_radiodialysis.R", "reactor_decision_tree.R", "scoby_3d.jl",
    "validate_serial.jl", "export_checkpoint.jl",
]
SIM_FILES = [joinpath(REPO, f) for f in SIM_FILE_NAMES]

"""
    _scan(paths, pattern) -> Vector{(path, lineno, line)}

Every line in every file matching `pattern` (a case-insensitive Regex). Reused by both the
real-file checks and their synthetic negative controls, so a negative control exercises the
same detection logic the real check trusts, not a reimplementation of it.
"""
function _scan(paths, pattern::Regex)
    hits = Tuple{String,Int,String}[]
    for p in paths
        for (i, line) in enumerate(eachline(p))
            occursin(pattern, line) && push!(hits, (p, i, line))
        end
    end
    hits
end

# Synthetic-file variant: same detection function, in-memory content instead of a path, so
# the negative controls below prove `_scan`'s REGEX actually matches something before any
# test trusts it finding nothing in the real files (AGENTS.md rule 1).
function _scan_text(name_content_pairs, pattern::Regex)
    hits = Tuple{String,Int,String}[]
    for (name, content) in name_content_pairs
        for (i, line) in enumerate(split(content, '\n'))
            occursin(pattern, line) && push!(hits, (name, i, line))
        end
    end
    hits
end

# ---------------------------------------------------------------- F_s identifier

let
    control = _scan_text([("fake.jl", "x = 1\nF_s = compute_fitness()\n")], r"\bF_s\b"i)
    @test length(control) == 1   # the detector can find a planted hit before it's trusted

    real = _scan(SIM_FILES, r"\bF_s\b"i)
    @test isempty(real)
end

# ---------------------------------------------------------------- "fitness" word, exact count

let
    # Both directions of negative control, since this is the one assertion that is an exact
    # count rather than an absence: proving the detector finds a hit doesn't prove it counts
    # correctly. A three-occurrence fixture must fail (over-count), and a one-occurrence
    # fixture must also fail (under-count) -- passing only one direction would leave a
    # detector that silently tolerates the other kind of drift.
    two_occurrences = [("a.jl", "# fitness\n"), ("b.R", "# fitness\n")]
    three_occurrences = [("a.jl", "# fitness\n"), ("b.R", "# fitness\n"), ("c.jl", "# fitness\n")]
    one_occurrence = [("a.jl", "# fitness\n"), ("b.R", "not here\n")]

    @test length(_scan_text(two_occurrences, r"fitness"i)) == 2
    @test length(_scan_text(three_occurrences, r"fitness"i)) != 2
    @test length(_scan_text(one_occurrence, r"fitness"i)) != 2

    real = _scan(SIM_FILES, r"fitness"i)
    @test length(real) == 2
    # every real match is a comment line (a header attribution, not an implemented quantity)
    @test all(occursin(r"^\s*#", line) for (_, _, line) in real)
    # and both live in exactly the two files the manuscript now names
    matched_files = Set(basename(p) for (p, _, _) in real)
    @test matched_files == Set(["biofilms_potts.jl", "biofilms_radiodialysis.R"])
end

# ---------------------------------------------------------------- momentum identifier
#
# NOT "momentum or velocity": biofilms_3d.R has an array its own comment calls a "velocity
# state" (reflection bookkeeping for the cylindrical wall boundary), but it is recomputed in
# full from the current position and forces at every step and never advanced by its own
# equation of motion -- it is not the persisted inertial degree of freedom the manuscript's
# claim is actually about. Checking for the bare word "velocity" would fail against code that
# does not violate the claim; the corrected .tex prose (section 3.4) makes this same
# distinction explicit. "momentum" has no such false-positive source and is checked directly.

let
    control = _scan_text([("fake.R", "momentum <- p * v\n")], r"\bmomentum\b"i)
    @test length(control) == 1

    real = _scan(SIM_FILES, r"\bmomentum\b"i)
    @test isempty(real)
end

# ---------------------------------------------------------------- unused bibliography entries

"""
    _unused_bibitem_keys(tex_text) -> Set{String}

Every `\\bibitem{key}` not reached by any `\\cite{...}` (comma-separated lists included).
"""
function _unused_bibitem_keys(tex_text::AbstractString)
    # STRIP LaTeX COMMENTS FIRST. A `\cite{foo}` inside a `%` comment counted as a
    # citation, so an orphan could be hidden behind one and never reported. The
    # Python side's normalise_markup has stripped `%` to end-of-line for exactly
    # this reason since it was written; this function did not.
    stripped = replace(tex_text, r"(?m)(?<!\\)%.*" => "")
    defined = [m.captures[1] for m in eachmatch(r"\\bibitem\{([^}]*)\}", tex_text)]
    cited = Set{String}()
    for m in eachmatch(r"\\cite\{([^}]*)\}", stripped)
        for key in split(m.captures[1], ',')
            push!(cited, strip(key))
        end
    end
    Set(k for k in defined if !(k in cited))
end

let
    # Detector correctness first, on a synthetic fixture with one used and one unused key --
    # @test_broken against the real manuscript below proves today's count matches what's
    # expected, not that the detector itself is right.
    fixture = "\\bibitem{used}\nSomething, cited.\n\\bibitem{unused}\nSomething, not cited.\n" *
              "Text with \\cite{used} in it.\n"
    @test _unused_bibitem_keys(fixture) == Set(["unused"])

    tex_path = joinpath(REPO, "preprint",
                         "modeling_radioresistance_and_radiotropic_fitness.tex")
    tex_text = read(tex_path, String)
    unused = _unused_bibitem_keys(tex_text)
    if !isempty(unused)
        println("  unused \\bibitem keys ($(length(unused))): ", join(sort(collect(unused)), ", "))
    end
    # KNOWN, PRE-EXISTING, NOT FIXED IN THIS PASS: 15 defined-but-unused entries, a separate
    # finding from PP-REF-* cleanup work, not something this change is responsible for.
    #
    # PIN THE COUNT, because @test_broken alone does not. It fails on zero and passes on
    # every other number, so 15 drifting to 16 -- a NEW unused entry, exactly what this
    # testset exists to prevent -- is indistinguishable from the status quo, and 15 dropping
    # to 14 loses the record of what was left. The ordinary assertion is the one with teeth
    # in both directions; @test_broken stays alongside it so reaching zero still reports an
    # UNEXPECTED PASS and forces someone to delete this whole block rather than let the
    # finding rot.
    # TWO CHANGES HERE FOR TWO DIFFERENT REASONS, AND THE COMMIT MUST NOT CONFLATE
    # THEM.
    #
    # (1) THE ASSERTION IS ON THE SET, NOT ON A COUNT, AND THAT IS A REPAIR
    # INDEPENDENT OF ANY NUMBER MOVING. `_unused_bibitem_keys` already RETURNS a
    # Set and the previous assertion threw it away to compare `length(...)` to a
    # constant. Cite one entry and orphan another in the same commit and the count
    # is unchanged and the test passes -- membership changed and nothing could see
    # it. That is the count-where-a-set-is-meant shape this repository has now had
    # caught three times, sitting in the guard for the bibliography.
    #
    # (2) The membership shrank from 15 to 13 because `turick2011` and
    # `casadevall2017` are now cited in section 2.6 -- as content, not as
    # bibliography tidying: the first is the one positive measurement in that
    # section and it runs against the mechanism, the second is the proponents
    # conceding the mechanism is undiscovered. A future shrink because somebody
    # deleted a bibitem to clear a red is a DIFFERENT event, and only the named
    # form below can tell the two apart.
    # THE SET CARRIES A CLASSIFICATION, NOT A LIST OF KEYS, so that a change in
    # membership says WHICH KIND of change it was. A bare key list makes "cited a
    # gap" and "deleted a bibitem to clear a red" the same edit. Resolved against
    # Crossref with PubMed/Europe PMC as second index on 2026-08-31; every entry
    # below had its DOI checked, not assumed, because four DOIs in this
    # bibliography have now resolved to unrelated papers or to nothing.
    #
    # (i) SUPPORTS AN UNCITED CLAIM -- a real gap, with the sentence identified.
    #     Citing these is a manuscript decision, not bibliography tidying, and two
    #     of them are corrections rather than footnotes. NOT applied here.
    # Four entries left this set on 2026-08-31 by being CITED, each for its own
    # reason and none of them bibliography tidying: alpkvist2007 and eberl2001
    # because the sentence that misattributed their content to xavier2005 was
    # corrected, and khajo2011 and malo2018 because section 2.6 cannot rest its
    # argument on naming what runs against it while omitting the two positives its
    # own audit files.
    # graner1992 LEFT THIS SET ON 2026-09-01 BY BEING CITED, on the Methods sentence
    # describing the copy-attempt construction. An origin cite for a named formalism is a
    # PROVENANCE claim -- the paper's title states the construction -- so it needed no
    # endpoint reading, which is what separated it from the two below.
    UNUSED_GAP = Set([
        # robertson2012's ENTRY IN THIS SET USED TO SAY "primary source for L228-229",
        # AND THAT WAS WRONG. Read from PubMed on 2026-09-01 (PMID 23139812), Robertson
        # 2012 reports that ionizing radiation "enhanced cell growth by increasing cell
        # division and cell size" and that low-dose exposure "significantly increased
        # survivability of BOTH the wild-type and the wdpks1 mutant" -- a positive growth
        # result, not a null -- with ribosomal biogenesis up in the irradiated wild type
        # but NOT in the mutant, a melanin-DEPENDENT difference. L228-229 claims comparable
        # NULL results across several modalities, and its actual primaries are the three
        # studies in radiotrophic_compatibility_audit.md 9.2 (GSE142318, GSE152116,
        # Microbiol Spectr 2023), reached through audit2026. Citing robertson2012 there
        # would have inserted a growth-enhancement paper as a null result: the section 2.6
        # defect, one key over, caught by reading the endpoint before writing the cite.
        # robertson2012 LEFT THIS SET ON 2026-09-01 BY BEING CITED IN 2.1, after the FULL TEXT
        # was read (PMC3490873) rather than the abstract. The abstract was not enough: it says
        # "we confirmed that ionizing radiation enhanced cell growth", and whether that was
        # measured or restated is the Khajo 2011 distinction. The Results heading settles it --
        # "Low Dose Ionizing Radiation Increases W. dermatitidis Growth Rate and Cell Size" --
        # and the same section reports the increase in the albino wdpks1 mutant too. So it is a
        # positive growth result whose melanin attribution fails, which is what 2.1 needed and
        # is NOT what L228-229 needed. Two readings, two different dispositions.
        # blasius1999 LEFT THIS SET ON 2026-09-01. Its gate was PP-T2-25's unapplied delete
        # of the Table 2 omega_s row; that delete is now applied, so the gate is gone and the
        # entry is cited in §2.3 where the Kuramoto lineage is discussed. Verified against
        # CROSSREF before citing, not reconstructed from the key: 10.1038/20676 returns
        # Blasius, Huppert & Stone (1999), Nature 399, 354-359, matching the bibitem on
        # title, all three authors, container, year, volume and pages. The key could equally
        # have been the Blasius boundary layer; reading it is what settled that.
    ])
    # (ii) DRAFT RESIDUE and (iii) DELIBERATE CONTEXT are both EMPTY as of 2026-09-01: the
    # six entries were DELETED from the bibliography, which is the disposition the
    # classification already recorded for them. brim2000 and kazy2009 supported claims that
    # are gone; battista1997, eisenman2012, lloyd2005 and newsome2014 supported no specific
    # sentence. The sets stay named rather than being removed, so that a future entry
    # arriving in either category lands somewhere that says what the category means.
    UNUSED_RESIDUE = Set(String[])
    UNUSED_CONTEXT = Set(String[])
    KNOWN_UNUSED = union(UNUSED_GAP, UNUSED_RESIDUE, UNUSED_CONTEXT)

    if unused != KNOWN_UNUSED
        for (name, s) in (("gap", UNUSED_GAP), ("residue", UNUSED_RESIDUE),
                          ("context", UNUSED_CONTEXT))
            left = setdiff(s, unused)
            !isempty(left) && println("  no longer unused, was $name: ",
                                      join(sort(collect(left)), ", "))
        end
        println("  newly unused: ",
                join(sort(collect(setdiff(unused, KNOWN_UNUSED))), ", "))
    end
    @test unused == KNOWN_UNUSED

    # THE ASSERTION MUST CONSUME THE FINEST-GRAINED THING THE FUNCTION RETURNS.
    # The line above compares the UNION, so moving graner1992 from UNUSED_GAP (a
    # real citation gap) to UNUSED_CONTEXT ("dropping it loses nothing citable")
    # -- opposite meanings -- leaves the union unchanged and the suite green. The
    # categories were printed on mismatch and never compared, which made the
    # classification documentation rather than an assertion, in the guard built
    # so that a decrement could not hide a category change.
    #
    # THIRD INSTANCE OF ONE SHAPE: count-where-a-set-was-meant (this same
    # assertion, earlier), subset-where-full-coverage-was-meant (the strict row
    # parser in test_guide_citations.py), and now
    # categories-printed-where-membership-was-meant. Every time the discriminating
    # data was already computed and discarded at the assertion line.
    @test isdisjoint(UNUSED_GAP, UNUSED_RESIDUE)
    @test isdisjoint(UNUSED_GAP, UNUSED_CONTEXT)
    @test isdisjoint(UNUSED_RESIDUE, UNUSED_CONTEXT)
    @test "graner1992" ∉ UNUSED_GAP        # cited 2026-09-01, so it must have LEFT the set
    @test "robertson2012" ∉ UNUSED_GAP     # cited 2026-09-01 in §2.1, after the FULL TEXT
    @test "blasius1999" ∉ UNUSED_GAP       # cited 2026-09-01, gate removed with PP-T2-25
    @test isempty(UNUSED_GAP) && isempty(UNUSED_RESIDUE) && isempty(UNUSED_CONTEXT)

    # FLIPPED FROM @test_broken ON 2026-09-01, AND THE TWO ARE NOT THE SAME CLAIM.
    # @test_broken said "we expect orphans and want to be told when there are none" -- it
    # tracked a goal. @test says "an uncited bibitem is a failure" -- it enforces a rule.
    # THE COST OF THE STRICT FORM IS REAL AND IS BEING ACCEPTED DELIBERATELY: the natural
    # workflow is add-the-bibitem-then-cite-it, and this makes that intermediate state red.
    # That is the intended trade -- every entry now in the bibliography is cited, and the
    # classification sets above exist to make a future exception a recorded decision rather
    # than a silent one. Someone who needs the intermediate state adds the key to the right
    # category set, which is a one-line edit that says which kind of orphan it is.
    @test isempty(unused)
end

# ------------------------------------------------- planned feedback is unimplemented
#
# Section 3.11 is the third "specification, not method" subsection, and like the other two
# its claim is checkable: nothing it describes exists in code. The capacity ceiling, the
# damage scalar and the bulk-water diffusivity limit are all absent -- but k_ads is NOT,
# because biofilms_radiodialysis.R really does carry a rate constant. The prose says "no
# adsorption capacity", meaning the ceiling, and this distinguishes the two rather than
# checking a looser word that would fail against correct code.

let
    control = _scan_text([("fake.R", "X_max <- 1.0\nq_max <- 2.0\n")], r"\b(X_max|q_max)\b"i)
    @test length(control) == 2   # the detector finds planted hits before it is trusted

    # The capacity ceiling and the damage scalar do not exist.
    for pattern in (r"\bX_max\b"i, r"\bq_max\b"i, r"\bLangmuir\b"i,
                    r"\bgamma_damage\b"i, r"\blambda_rad\b"i, r"\bD_w\b"i)
        @test isempty(_scan(SIM_FILES, pattern))
    end

    # ...but the RATE constant does, and the prose must not be read as denying it.
    # A test asserting "no k_ads anywhere" would fail against correct code, which is
    # how an over-broad absence claim gets discovered only after it ships.
    k_ads_hits = _scan(SIM_FILES, r"\bk_ads\b"i)
    @test !isempty(k_ads_hits)
    @test any(occursin("radiodialysis", basename(p)) for (p, _, _) in k_ads_hits)

    # And the sorption term in that file is reversible: k_des is present, so the
    # manuscript's claim that the implemented solver is not an irreversible sink is
    # checked rather than asserted.
    @test !isempty(_scan(SIM_FILES, r"\bk_des\b"i))
end

# ------------------------------------------- the exclusion, ASSERTED not inspected
#
# THE CHECK ABOVE IS SCOPED TO A CURATED NAME LIST, AND THAT CUTS BOTH WAYS.
# SIM_FILE_NAMES is nine hand-maintained entries, not a glob. So "X_max appears in
# no simulation file" stays green not only when no such file has one, but also when
# a file that has one simply is not on the list -- and nothing states which of those
# two is the case. That is the fixed-name-list hazard arriving from the far side: a
# guard passing for the wrong reason, which is the shape this file already records
# for the prose-broader-than-the-check defect one level up.
#
# It became live when analysis/henry_langmuir_bound.R was added. That producer
# necessarily contains X_max, q_max and Langmuir -- it exists to compare the forms
# -- and it is correctly outside SIM_FILE_NAMES because it is not a repo-root
# simulation source. But "correctly outside" was, until this block, a fact somebody
# confirmed once by reading, not a property anything enforced.
#
# So the enumeration is inverted: every file in the repository carrying any of the
# three terms must be on a DECLARED list. A new file gets added deliberately or the
# suite fails. The list is the claim; the scan is what checks it.
CEILING_VOCAB_ALLOWED = [
    # the manuscript states the forms and contrasts them -- section 3.12
    "preprint/modeling_radioresistance_and_radiotropic_fitness.tex",
    "preprint/wan_meeting_handout.tex",
    "calibration/tests/fixtures/wan_meeting_handout_prefix.tex",  # its known-bad
    # the register and the ledger record that the ceiling is absent
    "data/calibration/suspended_isotherm_proposal.csv",
    "data/calibration/sop_index.csv",
    "data/claims_ledger.csv",
    # the producer for PP-SORP-01, and this file's own assertions
    "analysis/henry_langmuir_bound.R",
    "tests/manuscript_claims_tests.jl",
]

function _ceiling_vocab_files(root)
    hits = String[]
    for (dir, _, files) in walkdir(root)
        occursin(joinpath(root, ".git"), dir) && continue
        for f in files
            p = joinpath(dir, f)
            rel = relpath(p, root)
            startswith(rel, ".git") && continue
            text = try
                read(p, String)
            catch
                continue        # binary or unreadable: cannot carry the vocabulary as text
            end
            if occursin(r"\bX_max\b"i, text) || occursin(r"\bq_max\b"i, text) ||
               occursin(r"\bLangmuir\b"i, text)
                push!(hits, replace(rel, '\\' => '/'))
            end
        end
    end
    sort(hits)
end

@testset "every file carrying capacity-ceiling vocabulary is declared" begin
    found = _ceiling_vocab_files(REPO)
    @test found == sort(CEILING_VOCAB_ALLOWED)

    # CONTROL: the enumeration must be able to find a file that is not declared.
    # Without this, an allow-list check passes identically when the walk is broken,
    # which is the always-passes half of a mutation harness wearing a set equality.
    planted = joinpath(REPO, "_ceiling_vocab_control.tmp")
    try
        write(planted, "q_max <- 1.0\n")
        @test "_ceiling_vocab_control.tmp" in _ceiling_vocab_files(REPO)
        @test _ceiling_vocab_files(REPO) != sort(CEILING_VOCAB_ALLOWED)
    finally
        rm(planted; force = true)
    end
    @test _ceiling_vocab_files(REPO) == sort(CEILING_VOCAB_ALLOWED)  # restored
end

# ---------------------------------------------- the acceptance path reads no nutrient
#
# §2.6 states that the CPM does not resolve carbon as a driver of the trajectory. That
# rests on a COUNT taken on one day -- zero nutrient reads inside two function bodies --
# and a counted zero expires exactly like the F_s grep above: nothing stops a later commit
# adding a nutrient term to the acceptance rule, at which point the manuscript's claim
# becomes false and the manuscript is the artifact that cannot notice.
#
# SCOPE, AS NARROW AS THE CODE ACTUALLY COUNTED. Not "the CPM ignores nutrients" and not
# "the field is unused" -- both are broader than what was checked and the second is FALSE.
# `state.nutrient` is initialised by `init_nutrient!`, integrated every step by
# `update_nutrient!`, serialised in export_checkpoint.jl, round-tripped in
# checkpoint_io_tests.jl and parity-checked in the JACC port. It is a live field that the
# ACCEPTANCE PATH does not consult, which is the claim, and calibration's
# spatial/time_observable.py records the same fact independently.

"Every column-0 `function NAME(...)` in `path`, as name => body lines."
function _function_bodies(path::AbstractString)
    lines = readlines(path)
    out = Dict{String,Vector{String}}()
    for (i, l) in enumerate(lines)
        m = match(r"^function\s+([A-Za-z_][A-Za-z0-9_!]*)\(", l)
        m === nothing && continue
        j = findnext(x -> x == "end", lines, i)
        j === nothing && continue
        out[m[1]] = lines[i:j]
    end
    out
end

"Body of `function NAME(` in `path`, up to the first column-0 `end`."
function _function_body(path::AbstractString, name::AbstractString)
    lines = readlines(path)
    i = findfirst(l -> startswith(l, "function $name("), lines)
    i === nothing && error("no `function $name(` at column 0 in $path")
    j = findnext(l -> l == "end", lines, i)
    j === nothing && error("unterminated `function $name` in $path")
    return lines[i:j]
end

let
    # A READ, not a mention: `state.nutrient` or `nutrient[...]`. A comment naming the
    # field must not fail this, or the check becomes a trap for the next person who
    # documents why the field is absent here.
    reads = r"\.nutrient\b|\bnutrient\s*\["i

    # THE DETECTOR FINDS A PLANTED READ BEFORE IT IS TRUSTED TO FIND NONE.
    planted = ["function f(x)", "    C = state.nutrient", "    y = nutrient[1, 2, 3]", "end"]
    @test count(l -> occursin(reads, l), planted) == 2
    # ...and a mention that is not a read does not trip it.
    @test !occursin(reads, "    # never reads the nutrient field")

    serial = joinpath(REPO, "biofilms_potts.jl")

    # THE SENTENCE SAYS "THE TRAJECTORY", SO THE CHECK COVERS EVERY FUNCTION THAT MOVES
    # THE LATTICE, not just the two whose bodies were first counted. Enforcing at one
    # scope while stating at another is the defect the F_s check above exists for.
    writes_lattice = r"^\s*(lat|lattice)\[[^\]]*\]\s*=[^=]"
    movers = String[]
    for (name, body) in _function_bodies(serial)
        any(l -> occursin(writes_lattice, l), body) && push!(movers, name)
    end

    # The detector must find a planted write before the set it produces is trusted.
    @test occursin(writes_lattice, "    lat[tx, ty, tz] = sigma")
    @test !occursin(writes_lattice, "    sigma = lat[tx, ty, tz]")   # a read is not a write

    # AND THE SET ITSELF IS PINNED. A new lattice mover appearing fails here, which forces
    # someone to check it rather than letting it inherit a claim made before it existed.
    # THIS LIST WAS WRONG WHEN WRITTEN BY HAND: it named divide_cell! and missed
    # place_cell!. divide_cell! mutates the lattice THROUGH place_cell! rather than
    # directly, so a manual grep for the writers produced a set that was wrong in both
    # directions, and this assertion is what said so on its first run.
    @test sort(movers) == sort(["compute_delta_H_terms", "init_state", "mcs_step!",
                                "place_cell!"])

    for fn in movers
        body = _function_body(serial, fn)
        @test !isempty(body)
        @test isempty(filter(l -> occursin(reads, l), body))
    end

    # THE FIELD IS LIVE, WHICH IS WHY THE CLAIM IS ABOUT THE ACCEPTANCE PATH AND NOT ABOUT
    # THE FIELD. If this ever came back empty the §2.6 sentence would be describing dead
    # state and would need rewriting, not re-passing.
    integrator = _function_body(serial, "update_nutrient!")
    @test any(l -> occursin(reads, l), integrator)
end

# ------------------------------------------------- the phase-locking thread agrees with itself
#
# SCOPE: PROSE paragraphs naming the phase-locking machinery -- H_kNN, Gamma_s, phi_s,
# "phase-lock" -- in the manuscript .tex. NOT the symbol table row, NOT the equations that
# define these symbols: a definition is not a claim that the model runs it, and widening this
# into the equation blocks would make every definition require a disclaimer. That exclusion is
# written here rather than left implicit, because "prose only" alone reads as a preference and
# the next person would widen it.
#
# SECTION 2.5 SAID THE FUNCTIONAL DID THE WORK; SECTIONS 3.2 AND 5 SAID IT WAS NEVER RUN, AND
# NOTHING COMPARED THEM. That was fixed in 567b135 -- and the guard written there matched
# `H_kNN` alone, so it did not catch §2.3 saying the Gamma_s kernel "draws directly on this
# tradition, extending it to multispecies microbial communities". Same defect, different
# symbol, invisible to a guard scoped to the symbol that prompted it. A GUARD WRITTEN FROM ONE
# INSTANCE MATCHES ONE INSTANCE; this one is scoped to the class.
#
# THE ASSERTION IS AGREEMENT, NOT ABSENCE OF A PHRASE. Banning a withdrawn verb is a guard the
# next synonym walks past. What must hold is that any paragraph claiming something about this
# machinery also says what its status is, in the vocabulary §3.2 and §5 already settled.
#
# PARAGRAPH-LEVEL, NOT SENTENCE-LEVEL, and that is load-bearing: §2.3 names Gamma_s in one
# sentence and carries "It remains unexercised" in the next, so a sentence-split check would
# fail on correct prose and invite someone to weaken it.
@testset "every prose paragraph naming the phase-locking machinery states its status" begin
    tex = read(joinpath(REPO, "preprint",
                        "modeling_radioresistance_and_radiotropic_fitness.tex"), String)

    names = r"Hamiltonian kNN|H_\{\\mathrm\{kNN\}\}|\\Gamma_s|\\phi_s|phase-lock"

    # MASK NON-PROSE FIRST, THEN MAP EACH PROSE HIT TO ITS SUBSECTION.
    # A definition is not a claim: the symbol-table row "$\Gamma_s$ & Phase-locked kernel" and
    # the equations that define these symbols must not require a disclaimer, or every
    # definition in the paper would. Masking preserves offsets so the containing subsection
    # can still be found.
    #
    # SUBSECTION GRANULARITY IS A CORRECTION MADE AFTER MEASURING. A first version checked
    # PARAGRAPHS and produced two false positives -- §3.4's leapfrog paragraph, whose
    # subsection carries "specification, not method" in a SIBLING paragraph, and §7.4's
    # structural-limits paragraph, which disclaims correctly with "absent from the executed
    # simulations". The manuscript marks status at subsection level, so a paragraph-level
    # guard cries wolf, and a guard that cries wolf gets weakened rather than obeyed.
    # The third thing it flagged was REAL and is why the widening happened: §2.2 said "Our
    # PSDE framework extends these precedents by incorporating ... phase-locking dynamics
    # that couple species interactions to external radiation fields" -- a fourth instance.
    masked = replace(tex,
        r"\\begin\{(equation|align|array|tabular)\}.*?\\end\{\1\}"s => (m -> " "^length(m)))
    masked = join([occursin(r"&.*\\\\", l) ? " "^length(l) : l
                   for l in split(masked, "\n")], "\n")

    heads = [m.offset for m in eachmatch(r"\\(sub)?section\{", tex)]
    owning(i) = maximum(h for h in heads if h <= i)
    function body(h)
        nxt = filter(x -> x > h, heads)
        tex[h:(isempty(nxt) ? lastindex(tex) : prevind(tex, minimum(nxt)))]
    end

    hits = unique(owning(m.offset) for m in eachmatch(names, masked))
    @test !isempty(hits)           # the detector must find something to be trusted
    @test length(hits) >= 4        # §2.2, §2.3, §3.2, §3.8, §5 at minimum

    # THE VOCABULARY IS DISCOVERED, NOT INVENTED: every phrase below was read out of a
    # subsection that already disclaims correctly. Adding one is a claim that the manuscript
    # uses it as a disclaimer, so it must be quoted from the manuscript.
    # ONE LITERAL, NOT THREE JOINED BY `*`. Regex `*` in Julia CONCATENATES: r"a|b" * r"c|d"
    # becomes (?:a|b)(?:c|d), which requires both in sequence. The first version of this line
    # was three alternation groups multiplied together and therefore matched almost nothing,
    # reporting four subsections as undisclaimed when only one was. It failed loudly, but it
    # failed as a FALSE POSITIVE -- the direction that gets a guard weakened.
    status = r"unexercised|unrepresentable|specified rather than implemented|specified but not|absent from the executed simulations|specification, not method|intended numerical treatment"
    for h in hits
        @test occursin(status, body(h))
    end

    @testset "the controls, one per symbol" begin
        # Each name must be detected, or a widened matcher that finds nothing passes silently.
        @test occursin(names, raw"Our phase-locking kernel $\Gamma_s(t,\mathbf{x})$ draws on")
        @test occursin(names, raw"the $H_{\mathrm{kNN}}$ machinery")
        @test occursin(names, raw"where $\phi_s(t) = \cos(\omega_s t - \theta_s)$")
        @test occursin(names, "Hamiltonian kNN Decision Tree")
        @test !occursin(names, "an unrelated paragraph about mesh convergence")

        # And the two withdrawn clauses must fail the status check, or it is decorative.
        @test !occursin(status, "a perspective that our Hamiltonian kNN decision tree " *
                                "operationalizes quantitatively.")
        @test !occursin(status, raw"Our phase-locking kernel $\Gamma_s(t,\mathbf{x})$ draws " *
                                "directly on this tradition, extending it to multispecies " *
                                "microbial communities.")
    end
end

# ------------------------------------------------- cross-references resolve, and are not numerals
#
# SCOPE: \ref and \label in the manuscript .tex, and literal `Section~<digit>` forms.
#
# FIGURE NUMERALS ARE DELIBERATELY NOT CHECKED HERE, AND THE REASON HAS AN OWNER.
# `Figure~1`, `Figures 1`, `Figure 2`, `Figure~3` and `Figure~4` are still literals in this
# file. They are claims_ledger row PP-FIG-01, verdict `delete`, required_to_fix "Delete or
# replace the references" -- UNAPPLIED. Acting on them is that row's decision and not this
# guard's, so widening this pattern to `Figure~` would silently execute a verdict nobody
# scheduled. Written here rather than in a commit message because a bare "scoped to Section~"
# reads as a design choice, and the next person reaching for this rule would widen it. That is
# how the __pycache__ exclusion went wrong: the reasoning was recorded somewhere the rule was
# not reached for.
#
# THE SECOND ASSERTION IS THE LOAD-BEARING ONE, BECAUSE CI CANNOT MAKE IT.
# manuscript-build runs `latexmk -pdf -interaction=nonstopmode -halt-on-error`. An undefined
# \ref is a LaTeX *warning*, not an error, so -halt-on-error does not stop, the job goes green,
# and the PDF ships the words "Section ??" -- a failure that looks like a typo to a reader and
# like success to the pipeline. Nothing else in the repository catches it. This file does.
@testset "section cross-references resolve, and none is a literal numeral" begin
    tex = read(joinpath(REPO, "preprint",
                        "modeling_radioresistance_and_radiotropic_fitness.tex"), String)

    numerals = [m.match for m in eachmatch(r"Section~[0-9]", tex)]
    @test isempty(numerals)

    labels = Set(m.captures[1] for m in eachmatch(r"\\label\{([^}]*)\}", tex))
    refs   = Set(m.captures[1] for m in eachmatch(r"\\ref\{([^}]*)\}", tex))
    dangling = sort(collect(setdiff(refs, labels)))
    @test isempty(dangling)          # a dangling \ref renders as "??" and CI stays green
    @test !isempty(refs)             # non-vacuity: the detector must find references at all

    @testset "the controls, in both directions" begin
        # Each assertion must fail on the form it was built from, or it is asserting over an
        # empty set and would pass on any manuscript whatsoever.
        @test !isempty([m.match for m in eachmatch(r"Section~[0-9]", "see Section~3.8, which")])
        @test isempty([m.match for m in eachmatch(r"Section~[0-9]",
                                                  raw"see Section~\ref{sec:knn}, which")])

        fake_labels = Set(["sec:intro"])
        fake_refs   = Set(["sec:intro", "sec:nonexistent"])
        @test !isempty(setdiff(fake_refs, fake_labels))       # a dangling ref is detected
        @test isempty(setdiff(Set(["sec:intro"]), fake_labels))
    end

    # The five converted references must each resolve, named so that deleting one is loud.
    for k in ("sec:intro", "sec:melanin_radiotrophy", "sec:knn", "sec:framework",
              "sec:discussion")
        @test k in labels
        @test k in refs
    end
end

# ------------------------------------------------- PP-RECIP-01's two sentences still exist
#
# SCOPE: the two manuscript sentences claims_ledger row PP-RECIP-01 connects. Nothing else.
#
# THIS GUARD IS NOT WHAT KEEPS THE ROW OFF THE ENFORCED NOWHERE REPORT, and saying so was
# wrong when it was proposed. That report's predicate is document-based --
# `r["document"] in PSEUDO_DOCUMENTS or GENERATED_ARTIFACTS`, and PSEUDO_DOCUMENTS is
# {"repository", "correspondence"} -- so a row whose document is `preprint` is off the report
# whatever its code_location says. (The corollary matters for the thirty-one rows that ARE on
# it: they cannot be cleared by adding guards, because they are listed for naming no file.)
#
# WHAT IT ACTUALLY BUYS: PP-RECIP-01 records that two separately-stated sentences are one
# mechanism. That is a claim ABOUT A RELATION, so it goes stale the moment either end is
# edited away -- and a record of a connection outliving the thing connected is the
# absence-record-going-stale family, caught here before it happens rather than after.
@testset "the two sentences PP-RECIP-01 connects are still in the manuscript" begin
    tex = read(joinpath(REPO, "preprint",
                        "modeling_radioresistance_and_radiotropic_fitness.tex"), String)

    mechanism = "radiolysis"                     # §3.11, the established descriptor
    reciprocity = "assumes dose--time reciprocity"   # §3.12, the caveat
    kinetics = "competing generation and recombination kinetics"

    @test occursin(mechanism, tex)
    @test occursin(reciprocity, tex)
    @test occursin(kinetics, tex)

    # The relation, not just the ends: §3.12's caveat must still be ABOUT radiolysis. If the
    # kinetics clause were rewritten without naming radiolysis products, PP-RECIP-01's claim
    # that these are one mechanism would be the thing that went stale.
    i = findfirst(reciprocity, tex)
    @test i !== nothing
    window = tex[first(i):min(lastindex(tex), first(i) + 600)]
    @test occursin("radiolysis", window)

    @testset "the control" begin
        # Each string must be absent from text that does not contain it, or these are
        # occursin calls that could not fail.
        @test !occursin(reciprocity, "an unrelated paragraph about mesh convergence")
        @test !occursin("radiolysis", "an unrelated paragraph about mesh convergence")
    end
end
