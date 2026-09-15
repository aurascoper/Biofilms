"""
Lexical census: the CPM sources must contain no executable reference to the independent
signal module or field identifiers.

SCOPE, because a guard's stated scope is a claim. What this certifies is exactly:
these five spellings, in the root-level executable sources `census_sources` walks,
outside `#` comments. It cannot see a renamed import (`import .InertSignal: update! as
bump`), an alias (`const SF = signal_field`), or entry through a struct field it does
not name. Section 6.7 claims something broader -- no signal quantity enters the CPM
state, Hamiltonian or RNG -- and a name list cannot establish that. The structural
form would assert the quantity set `compute_delta_H_terms` and `mcs_step!` read, which
fails on a new coupling whatever it is spelled; not built. Production also hashes the
immutable parent, which is the check that does not depend on spelling at all.

BLIND SPOTS, measured rather than supposed, and all five pinned in
tests/signal_field_tests.jl so a future edit to the stripping cannot move them silently:
  - a `#` inside a string literal truncates the line, hiding a real call after it
    (false negative, and the one that matters);
  - `#= ref =#` on a single line is missed (false negative);
  - a reference on an interior line of a multi-line `#= ... =#` block IS reported,
    since only `#` is split on -- a false positive, the opposite direction.

Every identifier the census names has a planted-file control in test_guard_identifiers.jl;
a guarded name with no control is a name nobody has watched refuse.
"""
function signal_references(path)
    isfile(path) || error("missing serial source")
    src = read(path,String)
    isempty(src) && error("empty serial source")
    refs = Tuple{Int,String}[]
    # `update_signal!` is split out of the \b-terminated group on purpose. A single
    # trailing \b after the alternation can never hold for an alternative ending in
    # `!`: `!` is a non-word character, and a call is always followed by `(`, a space
    # or end of line, all non-word, so the boundary fails. `update_signal!(state)`
    # passed this census silently until test_guard_identifiers.jl planted it.
    pattern = r"\b(?:InertSignal|autoinducer|signal_field)\b|\bupdate_signal!|\.signal\b"
    for (n,line) in enumerate(split(src,'\n'))
        code = split(line,'#';limit=2)[1]
        occursin(pattern,code) && push!(refs,(n,strip(code)))
    end
    refs
end

"""
The executable Julia sources at the repository root: what a run actually executes.
Directories are outside this walk on purpose and the reason is different for each --
`viewer/` and `diagnostics/inert_signal/` own the signal field by design, and `tests/`
names these identifiers in order to test for them. The boundary is "the simulation's
own sources", not "the files that happen to be clean".
"""
census_sources(root = normpath(joinpath(@__DIR__, "..", ".."))) =
    sort([joinpath(root, f) for f in readdir(root) if endswith(f, ".jl")])

if abspath(PROGRAM_FILE) == @__FILE__
    # A zero-byte source is skipped by the property that earns it, not by name:
    # signal_references refuses an empty file because an empty file and a clean file
    # grep identically. tests/signal_field_tests.jl asserts which sources are empty.
    paths = isempty(ARGS) ? filter(p -> filesize(p) > 0, census_sources()) : ARGS
    for path in paths
        refs = signal_references(path)
        isempty(refs) || error("CPM source references inert signal: $path $refs")
    end
    println("No signal identifiers in $(length(paths)) root CPM sources (lexical scope).")
end
