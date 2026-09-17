#!/usr/bin/env julia
# Per-parcel melanin over a seed range, and a comparison against the Odin port's file.
#
#   julia diagnostics/melanin_ensemble/parcel_sweep.jl <out.csv> \
#        [--seeds 42:297] [--n 40] [--parcels 6] [--at 100]
#   julia diagnostics/melanin_ensemble/parcel_sweep.jl --compare <ours.csv> <theirs.csv>
#
# WHY THIS IS A SECOND FILE AND NOT A COLUMN IN sweep.jl. The sweep's committed CSV and
# the receipts beside it are pinned by their own hashes. A new column would rewrite those
# files without changing one number in them. This writes its own file instead.
#
# WHY IT RUNS THE MODEL AGAIN RATHER THAN READING A TRAJECTORY. A Snapshot keeps only
# per-species totals, so no per-cell mean can be recovered from one. The state at the read
# MCS is the only thing with the parcels still in it, and `run_simulation` returns the
# final state. A run to exactly `--at` therefore yields both in one call.
#
# THE RNG STREAM IS THE SAME ONE. `run_simulation` seeds MersenneTwister(seed), calls
# init_state(params; seed = seed), then mcs_step! once per sweep. A run to 100 visits the
# same states as the first 100 steps of a run to 400. The proof is in this file's own
# output: its mel_site column must equal the committed sweep's mean_melanin, and the
# per-seed control below refuses if it does not.

using Printf

include(joinpath(@__DIR__, "parcel_melanin.jl"))

"""
Load the serial model, minus its figure-export half.

Copied in shape from `sweep.jl`, not included from it: including that file would define a
second `main` in this one's namespace, and a reader would have to know which one runs.
"""
function load_serial()
    src = read(joinpath(@__DIR__, "..", "..", "biofilms_potts.jl"), String)
    parts = split(src, "#  13. Figure export")
    length(parts) == 2 || error("expected exactly one figure-export split marker, found $(length(parts)-1)")
    M = Module(:SerialRef)
    Base.eval(M, :(using LinearAlgebra, Statistics, Random, Printf))
    Base.include_string(M, parts[1], "biofilms_potts.jl")
    M
end

function getopt(args, flag, default)
    i = findfirst(==(flag), args)
    isnothing(i) ? default : args[i + 1]
end

function parse_seeds(spec::AbstractString)
    occursin(":", spec) || return parse.(Int, split(spec, ","))
    lo, hi = parse.(Int, split(spec, ":"))
    collect(lo:hi)
end

const COLUMNS = ("seed", "mcs", "species", "volume", "ncells", "mel_site", "mel_parcel")

"""
Read a parcel CSV into `(seed, mcs, species) => row`, finding columns by name.

Both sides of the comparison are read by this one function, so a column the port renamed
is refused rather than read positionally. The two melanin values are kept as the strings
the producer printed as well as as numbers: six-decimal string equality is the strictest
claim available, and a float tolerance would hide a producer that rounds differently.
"""
function read_parcels(path::AbstractString)
    rows = Dict{NTuple{3,Int},NamedTuple}()
    col = nothing
    for line in eachline(path)
        isempty(strip(line)) && continue
        startswith(line, "#") && continue
        f = split(line, ",")
        if isnothing(col)
            idx = Dict(strip(String(f[i])) => i for i in eachindex(f))
            absent = [c for c in COLUMNS if !haskey(idx, c)]
            isempty(absent) || error("$path: header lacks $(join(absent, ", ")): $line")
            col = Dict(c => idx[c] for c in COLUMNS)
            continue
        end
        key = (parse(Int, f[col["seed"]]), parse(Int, f[col["mcs"]]), parse(Int, f[col["species"]]))
        haskey(rows, key) && error("$path: duplicate row for seed/mcs/species $key")
        rows[key] = (volume = parse(Int, f[col["volume"]]),
                     ncells = parse(Int, f[col["ncells"]]),
                     site_s = String(strip(f[col["mel_site"]])),
                     parcel_s = String(strip(f[col["mel_parcel"]])),
                     site = parse(Float64, f[col["mel_site"]]),
                     parcel = parse(Float64, f[col["mel_parcel"]]))
    end
    isnothing(col) && error("$path: no header line")
    isempty(rows) && error("$path: header but no data rows")
    return rows
end

"Compare two parcel CSVs and print the result. Returns true when every shared row agrees."
function compare(ours_path::AbstractString, theirs_path::AbstractString)
    ours = read_parcels(ours_path)
    theirs = read_parcels(theirs_path)
    shared = sort(collect(intersect(keys(ours), keys(theirs))))
    only_ours = length(ours) - length(shared)
    only_theirs = length(theirs) - length(shared)

    @printf("ours   %s: %d rows\n", ours_path, length(ours))
    @printf("theirs %s: %d rows\n", theirs_path, length(theirs))
    @printf("shared keys: %d   ours only: %d   theirs only: %d\n", length(shared), only_ours, only_theirs)
    isempty(shared) && error("no (seed, mcs, species) key appears in both files")

    agree = Dict(c => 0 for c in ("volume", "ncells", "mel_site", "mel_parcel"))
    worst = Dict("mel_site" => 0.0, "mel_parcel" => 0.0)
    examples = String[]
    for k in shared
        a, b = ours[k], theirs[k]
        a.volume == b.volume && (agree["volume"] += 1)
        a.ncells == b.ncells && (agree["ncells"] += 1)
        a.site_s == b.site_s && (agree["mel_site"] += 1)
        a.parcel_s == b.parcel_s && (agree["mel_parcel"] += 1)
        worst["mel_site"] = max(worst["mel_site"], abs(a.site - b.site))
        worst["mel_parcel"] = max(worst["mel_parcel"], abs(a.parcel - b.parcel))
        if a.parcel_s != b.parcel_s && length(examples) < 5
            push!(examples, @sprintf("  seed %d mcs %d species %d: ours %s theirs %s",
                                     k[1], k[2], k[3], a.parcel_s, b.parcel_s))
        end
    end

    n = length(shared)
    for c in ("volume", "ncells", "mel_site", "mel_parcel")
        @printf("%-11s %d of %d agree%s\n", c, agree[c], n,
                haskey(worst, c) ? @sprintf("   worst |diff| %.3e", worst[c]) : "")
    end
    isempty(examples) || (println("first disagreements on mel_parcel:"); foreach(println, examples))
    return agree["mel_parcel"] == n && agree["mel_site"] == n &&
           agree["volume"] == n && agree["ncells"] == n
end

function sweep_main(args)
    out = args[1]
    ispath(out) && error("destination already exists: $out")
    seeds   = parse_seeds(getopt(args, "--seeds", "42:297"))
    N       = parse(Int, getopt(args, "--n", "40"))
    parcels = parse(Int, getopt(args, "--parcels", "6"))
    at      = parse(Int, getopt(args, "--at", "100"))
    # Refused before the model loads and before the destination opens, so a bad request
    # leaves nothing behind. Same order of checks as sweep.jl, for the same reason.
    isempty(seeds) && error("no seeds in --seeds $(getopt(args, "--seeds", "42:297"))")
    allunique(seeds) || error("duplicate seed in --seeds $(getopt(args, "--seeds", "42:297"))")
    at > 0 || error("--at must be positive, got $at")

    SR = load_serial()
    # snapshot_interval = at, and the run stops at `at`, so the model takes a snapshot
    # there by both of its two conditions.
    p = Base.invokelatest(SR.CPMParams; N = N, n_cells_per_species = parcels,
                          snapshot_interval = at)
    names = SR.SPECIES_NAMES

    open(out, "w") do io
        println(io, "# per-parcel melanin sweep")
        println(io, "# N=$N parcels_per_species=$parcels read_at_mcs=$at")
        println(io, "# seeds=", join(seeds, " "))
        println(io, "# mel_site   = mean melanin over OCCUPIED SITES, volume-weighted. This is ",
                    "biofilms_potts.jl take_snapshot's mean_melanin, recomputed and checked against it.")
        println(io, "# mel_parcel = mean over PARCELS of each parcel's own mean, every parcel ",
                    "weighted equally. Definition supplied by F. Fink for the Odin port, 2026-09-16.")
        println(io, "# The two are different estimators of different quantities, not two ",
                    "computations of one.")
        println(io, "seed,mcs,species,volume,ncells,mel_site,mel_parcel")
        for seed in seeds
            t0 = time()
            state, traj = redirect_stdout(devnull) do
                Base.invokelatest(SR.run_simulation, p, at; seed = seed)
            end
            hits = findall(s -> s.mcs == at, traj)
            length(hits) == 1 ||
                error("seed $seed: expected exactly one snapshot at MCS $at, found $(length(hits))")
            snap = traj[hits[1]]
            sd_of = Dict(sd.species => sd for sd in snap.species_data)

            species_of = Dict(id => c.species for (id, c) in state.cells)
            r = parcel_means(state.lattice, state.melanin, species_of, SR.N_SPECIES)

            # THE CONTROL, RUN ON EVERY SEED. Reading a state is not the same act as the
            # model's own accounting, so the two must be made to meet. Equality is exact
            # because the visit order is identical; a tolerance here would accept a
            # reordering that silently changes which float sum is reported.
            for s in 1:SR.N_SPECIES
                sd = sd_of[s]
                r.mel_site[s] == sd.mean_melanin || error(
                    "seed $seed species $s: recomputed mel_site $(r.mel_site[s]) against the " *
                    "model's $(sd.mean_melanin). The lattice read disagrees with take_snapshot.")
                r.n_sites[s] == sd.total_volume || error(
                    "seed $seed species $s: counted $(r.n_sites[s]) sites against the registry's " *
                    "$(sd.total_volume). Volume bookkeeping disagrees with the lattice.")
                r.n_parcels[s] == sd.n_cells || error(
                    "seed $seed species $s: counted $(r.n_parcels[s]) parcels with sites against " *
                    "the registry's $(sd.n_cells) cells.")
            end

            for s in 1:SR.N_SPECIES
                @printf(io, "%d,%d,%d,%d,%d,%.6f,%.6f\n", seed, at, s,
                        r.n_sites[s], r.n_parcels[s], r.mel_site[s], r.mel_parcel[s])
            end
            @printf(stderr, "seed %d  %.1f s  (%s mel_site %.6f mel_parcel %.6f)\n",
                    seed, time() - t0, names[3], r.mel_site[3], r.mel_parcel[3])
        end
    end
    println("wrote ", out)
end

function main(args)
    isempty(args) && error("usage: parcel_sweep.jl <out.csv> [--seeds 42:297] [--n 40] " *
                           "[--parcels 6] [--at 100]   |   --compare <ours.csv> <theirs.csv>")
    if args[1] == "--compare"
        length(args) == 3 || error("usage: parcel_sweep.jl --compare <ours.csv> <theirs.csv>")
        exit(compare(args[2], args[3]) ? 0 : 1)
    end
    sweep_main(args)
end

# Parenthesised: `@__FILE__ && x` parses as `@__FILE__(&& x)`.
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
