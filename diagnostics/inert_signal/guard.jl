"""
Narrow static guard: the serial CPM source must contain no executable reference to
the independent signal module/field identifiers. This is a lexical census, not
a proof about arbitrary aliasing. Every identifier it names has a planted-file
control in test_guard_identifiers.jl; a guarded name with no control is a name
nobody has watched refuse. Production also hashes the immutable parent.
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

if abspath(PROGRAM_FILE) == @__FILE__
    path = isempty(ARGS) ? joinpath(@__DIR__,"..","..","biofilms_potts.jl") : only(ARGS)
    refs=signal_references(path)
    isempty(refs) || error("serial CPM references inert signal: $refs")
    println("No signal identifiers in serial executable source (lexical scope).")
end
