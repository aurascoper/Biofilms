#!/usr/bin/env julia
# Per-species mean radius over a seed range, under four arms of the acceptance rule.
#
#   julia diagnostics/melanin_ensemble/radial_sweep.jl <out.csv> \
#        [--seeds 42:297] [--n 40] [--parcels 6] [--at 100] \
#        [--arms full,no_rad,no_mel,null] [--reference sweep.csv]
#
# WHY THIS EXISTS. The manuscript's only radiotropism result is one run: seed 42, six
# parcels per species, 100 MCS, read as a mean radius per species. It says of itself that
# the exact null must ultimately be obtained by replicate simulations. This file is those
# replicates, and it runs the two
# radiation-derived terms of the Hamiltonian switched off one at a time so that the
# question "does either term move a parcel" is measured. The terms' magnitudes alone do not
# answer it.
#
# THE ARMS, in terms of the CPMParams fields they zero. Nothing in the model is edited.
#   full    the shipped parameters.
#   no_rad  beta_ion = 0 for every species: the direct radiation term dH_rad is 0 on
#           every move. The melanin path is untouched.
#   no_mel  alpha_M_species = 0 for every species: melanin is never produced, the field
#           stays at its initial 0, diffusion of 0 is 0, so dH_mel = 0.5 * M = 0 on every
#           move. This is the coupling switched off without touching the hard-coded 0.5
#           that tests/prose_bounds.jl reads out of compute_delta_H_terms: the trajectory
#           is the one a zero coupling would give, because dH_mel is the only place M
#           enters the acceptance rule.
#   null    both radiation-derived terms off. Adhesion and volume only. The MCS 0 row of
#           every arm is the SEEDING FRAME, a different thing: the parcels as placed,
#           before any move.
#
# THE ARMS SHARE THEIR INITIAL PLACEMENT. init_state seeds its own MersenneTwister(seed)
# for the parcel centres, and none of the zeroed fields is read there, so at MCS 0 every
# arm of a seed holds the same parcels at the same sites. The acceptance stream diverges
# after the first move whose dH differs. Comparisons across arms are therefore paired
# within a seed, as analyse.jl pairs species within a seed.
#
# WHAT mean_r IS. take_snapshot's mean_r: the mean over a species' parcels of each
# parcel's centre-of-mass distance from the z axis, every parcel weighted equally, in
# lattice units. mean_r_over_R divides by R = N/2, which is what Figure 1 plotted.
# This is a per-parcel mean; mean_melanin in the same row is a per-site mean. The
# README of this directory records that the two weight differently.
#
# WHY A NEW FILE. As parcel_sweep.jl says: the committed sweep CSVs and
# the receipts beside them are pinned by their own hashes.
#
# THE MACHINE IS IN THE HEADER. Julia's random streams match across platforms, but a
# floating-point result on Apple ARM and on x86 is not guaranteed bit-identical, and the
# manuscript will quote numbers from this file. julia_version and machine are written so
# that the file says where it came from.

using Printf

"""
Load the serial model, minus its figure-export half.

Copied in shape from sweep.jl, not included from it: including that file would define a
second `main` in this one's namespace.
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

# The arm vocabulary. An arm not in this table is refused before the model loads; a
# dispatch that defaulted to `full` for an unknown name would run the wrong experiment
# under the requested label.
const ARMS = ("full", "no_rad", "no_mel", "null")

"Which of the two radiation-derived vectors an arm zeroes."
zeroes(arm::AbstractString) = arm == "full"   ? (rad = false, mel = false) :
                              arm == "no_rad" ? (rad = true,  mel = false) :
                              arm == "no_mel" ? (rad = false, mel = true)  :
                              arm == "null"   ? (rad = true,  mel = true)  :
                              error("unknown arm $arm; the arms are $(join(ARMS, ", "))")

"Build the parameters for one arm. The shipped vectors are read from the model's own defaults."
function arm_params(SR, arm, N, parcels, at)
    z = zeroes(arm)
    base = Base.invokelatest(SR.CPMParams; N = N, n_cells_per_species = parcels,
                             snapshot_interval = at)
    nsp = SR.N_SPECIES
    β = z.rad ? zeros(nsp) : copy(Base.invokelatest(getfield, base, :β_ion))
    α = z.mel ? zeros(nsp) : copy(Base.invokelatest(getfield, base, :α_M_species))
    Base.invokelatest(SR.CPMParams; N = N, n_cells_per_species = parcels,
                      snapshot_interval = at, β_ion = β, α_M_species = α)
end

"""
Read a sweep.jl CSV into (seed, mcs, species) => mean_melanin, as the printed string.

Used for the reference check below. Only sweep.jl's columns are needed; they are found by
name, and the file's own header line says which configuration it was run at.
"""
function read_reference(path::AbstractString)
    rows = Dict{NTuple{3,Int},String}()
    meta = Dict{String,Int}()
    col = nothing
    for line in eachline(path)
        isempty(strip(line)) && continue
        if startswith(line, "#")
            for m in eachmatch(r"(N|parcels_per_species)=(\d+)", line)
                meta[m[1]] = parse(Int, m[2])
            end
            continue
        end
        f = split(line, ",")
        if isnothing(col)
            idx = Dict(strip(String(f[i])) => i for i in eachindex(f))
            for c in ("seed", "mcs", "species", "mean_melanin")
                haskey(idx, c) || error("$path: header lacks $c: $line")
            end
            col = idx
            continue
        end
        key = (parse(Int, f[col["seed"]]), parse(Int, f[col["mcs"]]), parse(Int, f[col["species"]]))
        haskey(rows, key) && error("$path: duplicate row for seed/mcs/species $key")
        rows[key] = String(strip(f[col["mean_melanin"]]))
    end
    isnothing(col) && error("$path: no header line")
    isempty(rows) && error("$path: header but no data rows")
    haskey(meta, "N") && haskey(meta, "parcels_per_species") ||
        error("$path: header does not state N and parcels_per_species")
    rows, meta
end

"Mean over a species' cells of the centre-of-mass radius, recomputed from the registry."
function registry_mean_r(SR, state, N)
    r_sum = zeros(Float64, SR.N_SPECIES)
    r_cnt = zeros(Int, SR.N_SPECIES)
    for (_, cell) in state.cells
        s = cell.species
        r_sum[s] += Base.invokelatest(SR.radial_dist, cell.com[1], cell.com[2], N)
        r_cnt[s] += 1
    end
    [r_cnt[s] > 0 ? r_sum[s] / r_cnt[s] : 0.0 for s in 1:SR.N_SPECIES]
end

function sweep_main(args)
    out = args[1]
    ispath(out) && error("destination already exists: $out")
    seeds   = parse_seeds(getopt(args, "--seeds", "42:297"))
    N       = parse(Int, getopt(args, "--n", "40"))
    parcels = parse(Int, getopt(args, "--parcels", "6"))
    at      = parse(Int, getopt(args, "--at", "100"))
    arms    = String.(split(getopt(args, "--arms", join(ARMS, ",")), ","))
    refpath = getopt(args, "--reference", nothing)
    # Refused before the model loads and before the destination opens, so a bad request
    # leaves nothing behind. Same order as parcel_sweep.jl, for the same reason.
    isempty(seeds) && error("no seeds in --seeds $(getopt(args, "--seeds", "42:297"))")
    allunique(seeds) || error("duplicate seed in --seeds $(getopt(args, "--seeds", "42:297"))")
    at > 0 || error("--at must be positive, got $at")
    isempty(arms) && error("no arms in --arms")
    allunique(arms) || error("duplicate arm in --arms $(join(arms, ","))")
    foreach(zeroes, arms)                     # an unknown arm errors here, by name
    # The reference is read before the model loads too: a missing or malformed file is a
    # bad request; it says nothing about seed 42.
    reference = nothing
    if !isnothing(refpath)
        rows, meta = read_reference(refpath)
        (meta["N"], meta["parcels_per_species"]) == (N, parcels) ||
            error("$refpath is N=$(meta["N"]), $(meta["parcels_per_species"]) parcels; " *
                  "this run is N=$N, $parcels parcels, so the reference cannot bind it")
        reference = rows
    end

    SR = load_serial()
    names = SR.SPECIES_NAMES
    R = N / 2.0
    params = Dict(arm => arm_params(SR, arm, N, parcels, at) for arm in arms)
    checked = 0

    open(out, "w") do io
        println(io, "# radial ensemble sweep")
        println(io, "# N=$N parcels_per_species=$parcels read_at_mcs=$at R=$R")
        println(io, "# arms=", join(arms, ","))
        println(io, "# seeds=", join(seeds, " "))
        println(io, "# julia_version=$(VERSION) machine=$(Sys.MACHINE)")
        println(io, "# full: shipped CPMParams. no_rad: beta_ion = 0 for every species. ",
                    "no_mel: alpha_M_species = 0 for every species, so the melanin field ",
                    "stays 0 and dH_mel is 0 on every move. null: both.")
        println(io, "# mean_r = take_snapshot's mean over a species' PARCELS of the ",
                    "centre-of-mass radius, every parcel weighted equally, lattice units. ",
                    "mean_r_over_R divides by R = N/2. mean_melanin is a per-SITE mean, ",
                    "as in sweep.jl; the two weight differently.")
        println(io, "# Every arm of a seed starts from the same placement: init_state draws ",
                    "the centres from its own MersenneTwister(seed) and reads neither vector.")
        println(io, "seed,arm,mcs,species,species_name,ncells,volume,mean_r,mean_r_over_R,mean_melanin")
        for seed in seeds
            t0 = time()
            placement = nothing
            for arm in arms
                p = params[arm]
                state, traj = redirect_stdout(devnull) do
                    Base.invokelatest(SR.run_simulation, p, at; seed = seed)
                end
                hits = findall(s -> s.mcs == at, traj)
                length(hits) == 1 ||
                    error("seed $seed arm $arm: expected exactly one snapshot at MCS $at, found $(length(hits))")
                first = findall(s -> s.mcs == 0, traj)
                length(first) == 1 ||
                    error("seed $seed arm $arm: expected exactly one snapshot at MCS 0, found $(length(first))")
                snap0, snap = traj[first[1]], traj[hits[1]]

                # THE CONTROLS, RUN ON EVERY SEED AND ARM.
                # 1. Reading the registry is a second act beside the model's own snapshot;
                #    the two are made to meet exactly (same visit: a Dict iterated once).
                rec = registry_mean_r(SR, state, N)
                for sd in snap.species_data
                    rec[sd.species] == sd.mean_r || error(
                        "seed $seed arm $arm species $(sd.species): recomputed mean_r $(rec[sd.species]) " *
                        "against the snapshot's $(sd.mean_r). The registry read disagrees with take_snapshot.")
                end
                # 2. An arm that zeroes production must report a melanin field of exactly 0.
                if zeroes(arm).mel
                    for sd in snap.species_data
                        sd.mean_melanin == 0.0 || error(
                            "seed $seed arm $arm species $(sd.species): mean_melanin $(sd.mean_melanin) " *
                            "with alpha_M = 0. The melanin field has a source this file does not know about.")
                    end
                end
                # 3. The arms share their placement: MCS 0 must agree across arms per seed.
                here = [(sd.species, sd.n_cells, sd.total_volume, sd.mean_r) for sd in snap0.species_data]
                if isnothing(placement)
                    placement = here
                else
                    here == placement || error(
                        "seed $seed arm $arm: the MCS 0 snapshot differs from the first arm's. " *
                        "The arms do not share a placement, and the pairing below would be false.")
                end
                # 4. Against the committed sweep, when one is given: the full arm's melanin
                #    at MCS `at` must print identically to sweep.jl's, which binds this
                #    producer to the stream the manuscript's other numbers came from.
                if arm == "full" && !isnothing(reference)
                    for sd in snap.species_data
                        key = (seed, at, sd.species)
                        haskey(reference, key) || continue
                        printed = @sprintf("%.6f", sd.mean_melanin)
                        printed == reference[key] || error(
                            "seed $seed species $(sd.species): full-arm mean_melanin $printed " *
                            "against the reference's $(reference[key]). This producer is not on " *
                            "the committed sweep's stream.")
                        checked += 1
                    end
                end

                for s in (snap0, snap), sd in s.species_data
                    @printf(io, "%d,%s,%d,%d,%s,%d,%d,%.6f,%.6f,%.6f\n", seed, arm, s.mcs,
                            sd.species, replace(names[sd.species], "," => ""),
                            sd.n_cells, sd.total_volume, sd.mean_r, sd.mean_r / R, sd.mean_melanin)
                end
            end
            @printf(stderr, "seed %d  %.1f s  (%d arms)\n", seed, time() - t0, length(arms))
        end
    end
    if !isnothing(reference)
        # A reference that matched nothing is a check that never ran, and would read as
        # "passed". The count is printed and a zero is refused.
        checked > 0 || error("--reference $refpath shares no (seed, MCS $at, species) row with this run; " *
                             "the stream check never ran")
        println("reference check: $checked rows of the full arm matched $refpath")
    end
    println("wrote ", out)
end

function main(args)
    isempty(args) && error("usage: radial_sweep.jl <out.csv> [--seeds 42:297] [--n 40] " *
                           "[--parcels 6] [--at 100] [--arms full,no_rad,no_mel,null] " *
                           "[--reference sweep.csv]")
    sweep_main(args)
end

# Parenthesised: `@__FILE__ && x` parses as `@__FILE__(&& x)`.
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
