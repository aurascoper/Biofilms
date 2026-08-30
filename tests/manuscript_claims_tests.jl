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
    defined = [m.captures[1] for m in eachmatch(r"\\bibitem\{([^}]*)\}", tex_text)]
    cited = Set{String}()
    for m in eachmatch(r"\\cite\{([^}]*)\}", tex_text)
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
    KNOWN_UNUSED = 15
    @test length(unused) == KNOWN_UNUSED
    @test_broken isempty(unused)
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
    for fn in ("compute_delta_H_terms", "mcs_step!")
        body = _function_body(serial, fn)
        @test !isempty(body)
        hits = filter(l -> occursin(reads, l), body)
        @test isempty(hits)
    end

    # THE FIELD IS LIVE, WHICH IS WHY THE CLAIM IS ABOUT THE ACCEPTANCE PATH AND NOT ABOUT
    # THE FIELD. If this ever came back empty the §2.6 sentence would be describing dead
    # state and would need rewriting, not re-passing.
    integrator = _function_body(serial, "update_nutrient!")
    @test any(l -> occursin(reads, l), integrator)
end
