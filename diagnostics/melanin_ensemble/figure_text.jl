# figure_text.jl -- the stdlib-only parts of figure.jl: the sweep reader and the two
# strings whose wording is a claim about the data. Kept apart from the CairoMakie code so
# the diagnostic suite, which runs with no project, can exercise them, and so the render
# cannot print a sentence the data does not support.

"""
    read_at(path, at) -> (rows, alpha, meta)

Rows of a sweep CSV at MCS `at`, as seed => (species => mean_melanin); the declared
alpha_M per species; and the `# N=... parcels_per_species=... n_mcs=...` line as a Dict.
A repeated (seed, species) measurement at that MCS is refused: a concatenated or
duplicated sweep would otherwise be plotted as whichever copy came last, and equal
coefficients do not make two copies one run.
"""
function read_at(path, at)
    rows = Dict{Int,Dict{Int,Float64}}()
    alpha = Dict{Int,Float64}()
    meta = Dict{String,Int}()
    col = nothing
    for line in eachline(path)
        if startswith(line, "# N=")
            # sweep.jl's configuration line: "# N=40 parcels_per_species=6 n_mcs=400 ..."
            for kv in split(line[3:end])
                k, v = split(kv, "=")
                meta[k] = parse(Int, v)
            end
            continue
        end
        startswith(line, "#") && continue
        f = split(line, ",")
        if isnothing(col)
            col = Dict(name => i for (i, name) in enumerate(f))
            all(haskey(col, c) for c in ("seed", "mcs", "species", "alpha_M", "mean_melanin")) ||
                error("$path: header lacks a column this figure reads: $line")
            continue
        end
        parse(Int, f[col["mcs"]]) == at || continue
        sp = parse(Int, f[col["species"]])
        seed = parse(Int, f[col["seed"]])
        r = get!(rows, seed, Dict{Int,Float64}())
        # Assigning over an existing entry kept the last copy and said nothing.
        haskey(r, sp) && error("$path: duplicate row for seed $seed, species $sp, MCS $at")
        r[sp] = parse(Float64, f[col["mean_melanin"]])
        a = parse(Float64, f[col["alpha_M"]])
        # Keeping only the last value seen let a mixed-configuration file through.
        get(alpha, sp, a) == a || error("$path: alpha_M for species $sp changes within the file, " *
                                        "$(alpha[sp]) then $a")
        alpha[sp] = a
    end
    isempty(rows) && error("no rows at MCS $at in $path")
    for k in ("N", "parcels_per_species", "n_mcs")
        haskey(meta, k) || error("$path: no '# N=... parcels_per_species=... n_mcs=...' line; cannot state provenance")
    end
    return rows, alpha, meta
end

"""
    ordering_title(ordered, nseeds) -> String

Panel A's title. "Every seed descends" is said only when every seed does; otherwise the
title states the computed count and nothing more.
"""
ordering_title(ordered::Integer, nseeds::Integer) =
    (ordered == nseeds ? "Every seed descends: " : "") *
    "$(ordered) of $(nseeds) display the α_M ordering"

"""
    seed_set_label(seeds) -> String

The seeds actually plotted, exactly: a contiguous run prints as `lo:hi`, anything else as
the runs and singletons it is made of, `42:45,50,57`. Printing `min:max` for a sparse
set claims every seed in between was plotted.
"""
function seed_set_label(seeds)
    s = sort(unique(collect(Int, seeds)))
    isempty(s) && error("no seeds to label")
    parts = String[]
    i = 1
    while i <= length(s)
        j = i
        while j < length(s) && s[j + 1] == s[j] + 1
            j += 1
        end
        push!(parts, j == i ? string(s[i]) : "$(s[i]):$(s[j])")
        i = j + 1
    end
    join(parts, ",")
end

ordinal(n) = n == 2 ? "second" : n == 3 ? "third" :
             string(n, n % 10 == 1 && n % 100 != 11 ? "st" : n % 10 == 2 && n % 100 != 12 ? "nd" :
                       n % 10 == 3 && n % 100 != 13 ? "rd" : "th")

"""
    rank_title(rank, nseeds) -> String

Panel B's title: where the published seed's CN − AN difference ranks among the seeds.
Rank one is "the smallest", not "the smallest smallest", which is what composing an
ordinal with the superlative printed for a fixture whose published seed ranked first.
"""
rank_title(rank::Integer, nseeds::Integer) =
    rank == 1 ? "The published seed is the smallest of $(nseeds)" :
                "The published seed is the $(ordinal(rank)) smallest of $(nseeds)"
