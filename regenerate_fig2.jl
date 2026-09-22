#!/usr/bin/env julia
# Regenerate fig2_melanin_accumulation ONLY.
#
#   julia --project=. regenerate_fig2.jl
#   julia --project=. regenerate_fig2.jl --self-test    # the refusal must bite
#
# WHY THIS EXISTS. The committed figure printed an ordering attribution inside the image,
# which SEED42-SCOPE-01 requalified on 2026-09-17: the ordering fails in 19 of 256 seeds,
# so the input does not fix it. The row holds the wording. A caption inside an image is the
# sink a prose audit cannot reach (FIG-01 to FIG-07 record that failure), so the
# generator's text changed and the artifact has to be rebuilt to match it.
#
# WHY NOT `julia --project=. biofilms_potts.jl`. Its entry point is main_coupled(), which
# cannot run while RADIODIALYSIS: BLOCKED stands. export_figures writes figs 1 and 2
# together, so this renders into a temporary directory and copies exactly one figure
# across, as regenerate_fig3.jl does.
#
# WHAT THE FIGURE PLOTS IS PINNED, EVERY POINT OF IT. The committed reference is
# diagnostics/melanin_ensemble/sweep_n40_p6_seed42_int20.csv, written by this repository's
# own sweep.jl at the figure's configuration and snapshot interval. Every plotted snapshot
# of every species must match it, and the species set must be complete: a truncated
# reference that happens to hold one row would otherwise validate one number and pass the
# rest. The 16-seed sweep is read as well, so the two committed artifacts must agree at
# MCS 100. The script REFUSES to copy on any mismatch, and --self-test proves the refusal
# bites by mutating a value taken from the reference itself.
using Printf, SHA

const HERE   = @__DIR__
const OUT    = joinpath(HERE, "preprint", "figures")
const NAME   = "fig2_melanin_accumulation"
const ENS    = joinpath(HERE, "diagnostics", "melanin_ensemble")
const REF    = joinpath(ENS, "sweep_n40_p6_seed42_int20.csv")
const SWEEP  = joinpath(ENS, "sweep_n40_p6_seeds42-57.csv")
const MCS    = 100
const STEP   = 20
const SEED   = 42

M = Module(:Fig2Regen)
Base.eval(M, :(using LinearAlgebra, Statistics, Random, Printf))
Base.include(M, joinpath(HERE, "biofilms_potts.jl"))

const CPMParams      = Base.eval(M, :CPMParams)
const SPECIES_NAMES  = Base.eval(M, :SPECIES_NAMES)
const run_simulation = Base.eval(M, :run_simulation)
const export_figures = Base.eval(M, :export_figures)

"""(mcs, species name) => mean melanin, for seed SEED, from a committed sweep CSV."""
function reference(path::String)
    out = Dict{Tuple{Int,String},Float64}()
    for line in eachline(path)
        (startswith(line, "#") || startswith(line, "seed,")) && continue
        f = split(line, ",")
        parse(Int, f[1]) == SEED || continue
        out[(parse(Int, f[2]), String(f[4]))] = parse(Float64, f[8])
    end
    return out
end

"""Every way the reference could be incomplete, named rather than assumed."""
function reference_gaps(ref::Dict{Tuple{Int,String},Float64})
    gaps = String[]
    for mcs in 0:STEP:MCS, name in SPECIES_NAMES
        haskey(ref, (mcs, name)) || push!(gaps, "MCS $mcs, $name")
    end
    extra = setdiff(Set(k[2] for k in keys(ref)), Set(SPECIES_NAMES))
    isempty(extra) || push!(gaps, "names the model does not have: $(sort(collect(extra)))")
    return gaps
end

"""Every plotted point that differs from the reference."""
function mismatches(traj, ref::Dict{Tuple{Int,String},Float64})
    bad = String[]
    plotted = [s for s in traj if s.mcs % STEP == 0 && s.mcs <= MCS]
    [s.mcs for s in plotted] == collect(0:STEP:MCS) ||
        push!(bad, "the plotted snapshots are $([s.mcs for s in plotted]), not $(collect(0:STEP:MCS))")
    for s in plotted, sd in s.species_data
        want = ref[(s.mcs, sd.name)]
        abs(sd.mean_melanin - want) <= 5e-7 ||
            push!(bad, @sprintf("MCS %d, %s: rendered %.6f, committed %.6f", s.mcs, sd.name, sd.mean_melanin, want))
    end
    return bad
end

function render()
    params = CPMParams(N = 40, n_cells_per_species = 6, snapshot_interval = STEP)
    _, traj = run_simulation(params, MCS; seed = SEED)
    return params, traj
end

"""Every reason the second sweep fails to confirm the first, named rather than assumed. A
missing row, a non-finite value and a disagreement are three different artifacts, and a
comparison written as `> tol` would read NaN as agreement, because NaN compares false."""
function overlap_problems(overlap::Dict{Tuple{Int,String},Float64}, ref::Dict{Tuple{Int,String},Float64})
    out = String[]
    for n in SPECIES_NAMES
        if !haskey(overlap, (MCS, n))
            push!(out, "no MCS $MCS row for $n")
        elseif !isfinite(overlap[(MCS, n)])
            push!(out, "$n has a non-finite value at MCS $MCS: $(overlap[(MCS, n)])")
        elseif !(abs(overlap[(MCS, n)] - ref[(MCS, n)]) <= 5e-7)
            push!(out, @sprintf("%s disagrees at MCS %d: %.6f against %.6f", n, MCS, overlap[(MCS, n)], ref[(MCS, n)]))
        end
    end
    return out
end

"""Render the four artifacts into `dir` and check them there. Nothing outside `dir` is
touched, so a failed extraction cannot leave a new image beside an old checksum."""
function build_figure(traj, params, dir::String; extract = default_extract)
    export_figures(traj, nothing, params; outdir = dir)
    pdf, txt = joinpath(dir, "$NAME.pdf"), joinpath(dir, "$NAME.txt")
    extract(pdf, txt)
    isfile(txt) || error("REFUSING TO COPY: the extraction wrote no text sidecar")
    body = read(txt, String)
    isempty(strip(body)) && error("REFUSING TO COPY: the text sidecar came out empty, and an empty " *
                                  "file greps the same as a clean one")
    occursin("SEED42-SCOPE-01", body) ||
        error("REFUSING TO COPY: the text sidecar does not name the ledger row the caption cites, " *
              "so the extraction is not of the figure this run drew")
    write(joinpath(dir, "$NAME.sha256"), bytes2hex(SHA.sha256(read(pdf))) * "\n")
    return ["$NAME.pdf", "$NAME.png", "$NAME.sha256", "$NAME.txt"]
end

"""The default extraction. The claims guard scans this .txt rather than the PDF, so it is part
of the operation and not a printed follow-up."""
function default_extract(pdf::String, txt::String)
    Sys.which("pdftotext") === nothing &&
        error("REFUSING TO COPY: pdftotext is needed to refresh the text sidecar")
    run(`pdftotext -layout $pdf $txt`)
end

"""The production boundary: confirm the second sweep, compare every plotted point, build and
check all four artifacts in a staging directory, and only then replace the committed ones.
The self-test drives THIS with every known-bad, so a refusal that stopped refusing, or one
that half-wrote, could not pass it."""
function validate_and_copy(traj, params, ref::Dict{Tuple{Int,String},Float64},
                           overlap::Dict{Tuple{Int,String},Float64}, outdir::String;
                           extract = default_extract)
    problems = overlap_problems(overlap, ref)
    isempty(problems) || error("REFUSING TO COPY: the second committed sweep does not confirm the " *
                               "first:\n  " * join(problems, "\n  "))
    bad = mismatches(traj, ref)
    isempty(bad) || error("REFUSING TO COPY: the render does not match the committed reference:\n  " *
                          join(bad, "\n  ") * "\nThe figure would plot numbers no committed run supports.")
    staging = mktempdir(; cleanup = false)
    names = try
        build_figure(traj, params, staging; extract = extract)
    catch err
        rm(staging; recursive = true)
        rethrow(err)
    end
    for n in names
        cp(joinpath(staging, n), joinpath(outdir, n); force = true)
    end
    rm(staging; recursive = true)
    return sort(readdir(outdir))
end

"""Run the boundary with one known-bad and require a refusal that wrote nothing."""
function must_refuse(label, traj, params, ref, overlap; extract = default_extract)
    dir = mktempdir(; cleanup = false)
    threw = try
        validate_and_copy(traj, params, ref, overlap, dir; extract = extract)
        false
    catch err
        occursin("REFUSING TO COPY", sprint(showerror, err)) || rethrow(err)
        true
    end
    threw || error("SELF-TEST: $label was copied")
    isempty(readdir(dir)) || error("SELF-TEST: $label still wrote $(readdir(dir))")
    rm(dir; recursive = true)
    println("  self-test: $label is refused at the copy boundary, and nothing is written")
end

ref = reference(REF)
gaps = reference_gaps(ref)
isempty(gaps) || error("REFUSING: the reference in $REF does not cover $(length(gaps)) of the " *
                       "plotted points: $(first(gaps, 5))")

if "--self-test" in ARGS
    # Every known-bad is taken from an artifact, never typed, and every one goes through
    # validate_and_copy, the boundary the real run uses. Testing the helpers alone would
    # leave a bypassed refusal green.
    params, traj = render()
    overlap = reference(SWEEP)

    victim = (STEP, SPECIES_NAMES[3])
    moved = copy(ref)
    moved[victim] = ref[victim] + 1e-3
    moved[victim] == ref[victim] && error("SELF-TEST: the mutation did not apply at $victim")
    must_refuse("a moved reference value", traj, params, moved, overlap)

    dropped = filter(kv -> kv.first != (MCS, SPECIES_NAMES[1]), overlap)
    length(dropped) == length(overlap) - 1 || error("SELF-TEST: the overlap row was not dropped")
    must_refuse("a dropped overlap row", traj, params, ref, dropped)

    nan = copy(overlap)
    nan[(MCS, SPECIES_NAMES[1])] = NaN
    isnan(nan[(MCS, SPECIES_NAMES[1])]) || error("SELF-TEST: the NaN did not apply")
    must_refuse("a non-finite overlap value", traj, params, ref, nan)
    must_refuse("an extraction that writes nothing", traj, params, ref, overlap;
                extract = (pdf, txt) -> write(txt, ""))
    # Nonempty text from the wrong figure: the empty case above is refused one line earlier,
    # so without this the identity check could be deleted and every outcome would be the same.
    must_refuse("an extraction of another figure", traj, params, ref, overlap;
                extract = (pdf, txt) -> write(txt, "Mean melanin field value at occupied sites\n"))

    allowed = mktempdir(; cleanup = false)
    validate_and_copy(traj, params, ref, overlap, allowed) ==
        ["$NAME.pdf", "$NAME.png", "$NAME.sha256", "$NAME.txt"] ||
        error("SELF-TEST: the committed artifacts did not copy the figure: $(readdir(allowed))")
    rm(allowed; recursive = true)
    println("  self-test: the committed artifacts write all four files, so the boundary refuses by evidence")
    exit(0)
end

params, traj = render()
for mcs in 0:STEP:MCS, n in SPECIES_NAMES[[1, 3, 5]]
    @printf("  MCS %3d  %-20s committed %.6f\n", mcs, n, ref[(mcs, n)])
end
validate_and_copy(traj, params, ref, reference(SWEEP), OUT)

@printf("\n  every plotted point matches the committed reference\n")
@printf("  wrote %s.{pdf,png,txt,sha256} into %s\n", NAME, OUT)
@printf("  sha256 %s\n", strip(read(joinpath(OUT, "$NAME.sha256"), String)))
