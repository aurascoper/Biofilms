"""
Narrow static guard: the serial CPM source must contain no executable reference to
the independent signal module/field identifiers. This is a lexical census, not
a proof about arbitrary aliasing. Production also hashes the immutable parent.
"""
function signal_references(path)
    isfile(path) || error("missing serial source")
    src = read(path,String)
    isempty(src) && error("empty serial source")
    refs = Tuple{Int,String}[]
    pattern = r"\b(?:InertSignal|update_signal!|autoinducer|signal_field)\b|\.signal\b"
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
