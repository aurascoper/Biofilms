#!/usr/bin/env julia
# The melanin ensemble at the published configuration, as one figure.
#
# WHAT THIS SHOWS AND WHY IT IS DRAWN PAIRED. The three producers share one
# lattice per run, so their melanin values are not independent draws: the
# within-run correlation is +0.299 for CS-CN and -0.287 for CN-AN, opposite in
# sign, which means one pooled standard deviation would overstate one pair and
# understate the other. Three independent box plots would be that error drawn.
# Panel A therefore connects each seed's three values with a line -- a line that
# descends throughout IS a seed displaying the alpha_M ordering, so "16 of 16"
# is read off the picture rather than asserted beside it.
#
# Panel B is the CN-AN paired difference, sorted, because that is the pair the
# N=20 sweep could not resolve (12 of 16, sign test p = 0.077) and the one the
# published figure separates by 0.03.
#
# THE PUBLISHED SEED IS MARKED because it is not typical: seed 42's CN-AN gap is
# the third smallest of sixteen and the ensemble mean is 9.4x it. A reader who
# takes the published figure as representative of the model is reading its
# weakest draw for this pair.
#
# Not a manuscript artifact. It lives in diagnostics/ and writes its own .sha256
# (of the PNG: cairo stamps a CreationDate into the PDF, so two renders of one
# drawing hash differently there and identically here) and .txt (pdftotext
# -layout of the PDF) sidecars, so it is already compliant if it ever moves to
# preprint/figures/, where the staleness guards would apply to it.
#
#   julia --project=diagnostics/melanin_ensemble diagnostics/melanin_ensemble/figure.jl [sweep.csv]
#
# The root Project.toml does not declare CairoMakie; this directory's does, with a
# Manifest pinning the version that rendered the committed PNG.
#
# A sweep at any other configuration (N, parcels, MCS, or the alpha_M coefficients
# that give the x-axis its order) is refused, and so is one missing a producer row
# or the published seed: the provenance line below is printed from constants, and
# the CSV's own header and columns must agree with them.

using CairoMakie, Printf, SHA, Statistics

const HERE   = dirname(@__DIR__) |> dirname
const CSVIN  = length(ARGS) >= 1 ? ARGS[1] :
               joinpath(@__DIR__, "sweep_n40_p6_seeds42-57.csv")
const OUTBASE = joinpath(@__DIR__, "melanin_ensemble_n40_p6")

# Provenance, as constants rather than prose, so the line printed into the image
# cannot drift from what was plotted.
const N          = 40
const PARCELS    = 6
const N_MCS      = 400
const AT_MCS     = 100
const PUBLISHED  = 42          # the seed tests/fixtures/serial_seed42.csv pins
const PRODUCERS  = [(3, "C. sphaerospermum", 0.140),
                    (1, "C. neoformans",     0.100),
                    (5, "A. niger",          0.065)]

# Validated with the dataviz palette validator, light surface: worst adjacent
# CVD dE 24.7 (protan), normal-vision dE 33.6, both well above the >= 8 target.
const ENSEMBLE = colorant"#2a78d6"
const MARKED   = colorant"#eb6834"
const INK      = colorant"#2f2f2e"
const MUTED    = colorant"#6b6b68"

ordinal(n) = n == 1 ? "smallest" : n == 2 ? "second" : n == 3 ? "third" :
             string(n, n % 10 == 1 && n % 100 != 11 ? "st" : n % 10 == 2 && n % 100 != 12 ? "nd" :
                       n % 10 == 3 && n % 100 != 13 ? "rd" : "th")

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
        get!(rows, parse(Int, f[col["seed"]]), Dict{Int,Float64}())[sp] =
            parse(Float64, f[col["mean_melanin"]])
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

function main()
    # Refused before any output is written: nothing here may leave a fresh image beside
    # a stale sidecar, and the .txt needs pdftotext.
    isnothing(Sys.which("pdftotext")) &&
        error("pdftotext not found: the .txt sidecar could not be written, refusing to render")
    data, alpha, meta = read_at(CSVIN, AT_MCS)
    (meta["N"], meta["parcels_per_species"], meta["n_mcs"]) == (N, PARCELS, N_MCS) ||
        error("$CSVIN is N=$(meta["N"]), $(meta["parcels_per_species"]) parcels, $(meta["n_mcs"]) MCS; " *
              "this figure is drawn for N=$N, $PARCELS parcels, $N_MCS MCS, and its provenance line would lie")
    # The x-axis order IS the alpha_M order; a sweep run with other coefficients would be
    # drawn under a label that no longer describes it.
    for (id, name, a) in PRODUCERS
        isapprox(get(alpha, id, NaN), a; atol = 1e-6) ||
            error("$CSVIN: alpha_M for $name (species $id) is $(get(alpha, id, "absent")), " *
                  "this figure is drawn for $a")
    end
    haskey(data, PUBLISHED) || error("$CSVIN has no seed $PUBLISHED, the marked one")
    for (s, v) in data, (id, name, _) in PRODUCERS
        haskey(v, id) || error("$CSVIN: seed $s has no row for $name (species $id) at MCS $AT_MCS")
    end
    seeds = sort(collect(keys(data)))
    xs    = 1:length(PRODUCERS)
    gapc  = [(s, data[s][1] - data[s][5]) for s in seeds]      # CN - AN
    sort!(gapc, by = last)
    pubrank = findfirst(t -> t[1] == PUBLISHED, gapc)          # computed, not asserted
    ordered = count(s -> data[s][3] > data[s][1] > data[s][5], seeds)

    fig = Figure(size = (1000, 460), backgroundcolor = colorant"#fcfcfb")

    axA = Axis(fig[1, 1], xticks = (collect(xs), [p[2] for p in PRODUCERS]),
               ylabel = "mean melanin over occupied sites",
               title = "Every seed descends: $(ordered) of $(length(seeds)) display the α_M ordering",
               titlealign = :left, xgridvisible = false,
               ygridcolor = (:black, 0.06), leftspinevisible = false,
               topspinevisible = false, rightspinevisible = false,
               xticklabelrotation = 0.0, xticklabelsize = 11)
    for s in seeds
        s == PUBLISHED && continue
        lines!(axA, xs, [data[s][p[1]] for p in PRODUCERS];
               color = (ENSEMBLE, 0.45), linewidth = 2)
        scatter!(axA, xs, [data[s][p[1]] for p in PRODUCERS];
                 color = (ENSEMBLE, 0.55), markersize = 8)
    end
    lines!(axA, xs, [data[PUBLISHED][p[1]] for p in PRODUCERS];
           color = MARKED, linewidth = 3.5)
    scatter!(axA, xs, [data[PUBLISHED][p[1]] for p in PRODUCERS];
             color = MARKED, markersize = 12, strokecolor = colorant"#fcfcfb",
             strokewidth = 2)

    axB = Axis(fig[1, 2], ylabel = "C. neoformans − A. niger, paired within seed",
               xlabel = "seeds, sorted by that difference",
               title = "The published seed is the $(ordinal(pubrank)) smallest of $(length(seeds))",
               titlealign = :left, xgridvisible = false,
               ygridcolor = (:black, 0.06), leftspinevisible = false,
               topspinevisible = false, rightspinevisible = false,
               xticksvisible = false, xticklabelsvisible = false)
    hlines!(axB, [0.0]; color = (:black, 0.35), linewidth = 1.5)
    for (i, (s, g)) in enumerate(gapc)
        c = s == PUBLISHED ? MARKED : ENSEMBLE
        lines!(axB, [i, i], [0.0, g]; color = (c, 0.5), linewidth = 2)
        scatter!(axB, [i], [g]; color = c,
                 markersize = s == PUBLISHED ? 13 : 9,
                 strokecolor = colorant"#fcfcfb", strokewidth = s == PUBLISHED ? 2 : 0)
    end
    mg = mean(last.(gapc))
    hlines!(axB, [mg]; color = (ENSEMBLE, 0.8), linewidth = 2, linestyle = :dash)
    text!(axB, 0.5, mg; text = @sprintf(" ensemble mean %+.3f", mg),
          align = (:left, :bottom), color = ENSEMBLE, fontsize = 11)
    pubgap = data[PUBLISHED][1] - data[PUBLISHED][5]
    text!(axB, pubrank + 0.4, pubgap;
          text = @sprintf(" seed %d: %+.4f", PUBLISHED, pubgap),
          align = (:left, :bottom), color = MARKED, fontsize = 11)

    Legend(fig[2, 1],
           [LineElement(color = (ENSEMBLE, 0.55), linewidth = 2),
            LineElement(color = MARKED, linewidth = 3.5)],
           ["the other $(length(seeds) - 1) seeds, one line each",
            "seed $(PUBLISHED), the run the published figure and the golden fixture use"];
           orientation = :horizontal, framevisible = false, labelsize = 11,
           labelcolor = INK, tellheight = true)

    # Two lines, because one ran off the page -- and a provenance line that is
    # clipped is worse than none: the .txt sidecar would carry the truncation.
    Label(fig[3, 1:2],
          @sprintf("biofilms_potts.jl run_simulation via diagnostics/melanin_ensemble/sweep.jl  |  N=%d, %d parcels/species, %d MCS, seeds %d:%d, read at MCS %d\nobservable: volume-weighted mean melanin over occupied sites, NOT the mean of per-parcel means  |  α_M is a declared input, so an ordering displays it and does not measure it",
                   N, PARCELS, N_MCS, minimum(seeds), maximum(seeds), AT_MCS);
          fontsize = 9, color = MUTED, halign = :left, justification = :left,
          tellwidth = false)

    save(OUTBASE * ".pdf", fig)
    save(OUTBASE * ".png", fig; px_per_unit = 2)
    # The sidecars the header promises. Written here, by the same run, so a changed CSV
    # cannot leave a stale extraction and checksum beside a fresh image.
    write(OUTBASE * ".sha256", bytes2hex(sha256(read(OUTBASE * ".png"))) * "\n")
    run(pipeline(`pdftotext -layout $(OUTBASE * ".pdf") -`; stdout = OUTBASE * ".txt"))
    @printf("wrote %s.{pdf,png}\n  %d of %d seeds ordered; CN-AN mean %+.4f, seed %d %+.4f (rank %d)\n",
            OUTBASE, ordered, length(seeds), mg, PUBLISHED, pubgap, pubrank)
end

main()
