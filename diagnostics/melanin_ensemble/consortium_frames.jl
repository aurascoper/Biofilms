#!/usr/bin/env julia
# Occupied-site frames of one seed under two arms, for the consortium figure and video.
#
#   julia diagnostics/melanin_ensemble/consortium_frames.jl <outdir> [--seed 42] [--n 40]
#        [--parcels 6] [--mcs 100] [--arms full,null]
#
# WHY THIS EXISTS. The manuscript had no picture of the lattice, and the animations in the
# tree could not be used: the consortium GIF records no seed and is N = 20 with two parcels,
# and the 2024 files predate the radiation-sign fix and carry no run record. This writes the
# frames the ensemble configuration produces under one named seed, with the provenance the
# figure guards ask for, so that the picture is regenerable from the file it names.
#
# WHAT IT WRITES. One CSV per (arm, mcs) holding every occupied site as x,y,z,species,cell
# (gzip them before committing: 202 frames are 13 MB plain and about 2 MB compressed, and
# tools/render_consortium.py reads the .csv.gz),
# and a manifest naming the seed, the configuration, the arms as CPMParams fields zeroed, the
# git commit of the tree it ran in, the Julia version and the machine. The MCS 0 frame of
# every arm is identical: the arms share their placement (radial_sweep.jl checks this on
# every seed, and this file checks it on the one it writes).
#
# WHY NOT export_checkpoint.jl. That chain writes HDF5 through the root project, which the
# machine that produced the ensemble CSV could not instantiate; this loads the model through
# the same split-marker sandbox as radial_sweep.jl and writes plain text, so the frames can be
# reproduced with Julia alone, and rendered with Python alone.
#
# THE RNG STREAM IS THE ONE THE ENSEMBLE USED: run_simulation seeds MersenneTwister(seed) and
# init_state(params; seed), and the frames are read from the trajectory of one run to `--mcs`
# with snapshot_interval 1, so frame k is the state after k steps of that run.

using Printf

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

const ARMS = ("full", "no_rad", "no_mel", "null")
zeroes(arm) = arm == "full"   ? (rad = false, mel = false) :
              arm == "no_rad" ? (rad = true,  mel = false) :
              arm == "no_mel" ? (rad = false, mel = true)  :
              arm == "null"   ? (rad = true,  mel = true)  :
              error("unknown arm $arm; the arms are $(join(ARMS, ", "))")

function arm_params(SR, arm, N, parcels)
    z = zeroes(arm)
    base = Base.invokelatest(SR.CPMParams; N = N, n_cells_per_species = parcels, snapshot_interval = 1)
    nsp = SR.N_SPECIES
    β = z.rad ? zeros(nsp) : copy(Base.invokelatest(getfield, base, :β_ion))
    α = z.mel ? zeros(nsp) : copy(Base.invokelatest(getfield, base, :α_M_species))
    Base.invokelatest(SR.CPMParams; N = N, n_cells_per_species = parcels, snapshot_interval = 1,
                      β_ion = β, α_M_species = α)
end

"The tree's commit, with -dirty when the working tree differs, as the run artifacts record it."
function git_describe()
    try
        sha = strip(read(`git rev-parse --short HEAD`, String))
        dirty = !isempty(strip(read(`git status --porcelain`, String)))
        sha * (dirty ? "-dirty" : "")
    catch
        "unknown"
    end
end

"Write every occupied site of the lattice as x,y,z,species,cell."
function write_frame(path, SR, state)
    lat = state.lattice
    N = size(lat, 1)
    n = 0
    open(path, "w") do io
        println(io, "x,y,z,species,cell")
        @inbounds for z in 1:N, y in 1:N, x in 1:N
            σ = lat[x, y, z]
            σ > 0 || continue
            haskey(state.cells, Int(σ)) || continue
            @printf(io, "%d,%d,%d,%d,%d\n", x, y, z, state.cells[Int(σ)].species, σ)
            n += 1
        end
    end
    n
end

function main(args)
    isempty(args) && error("usage: consortium_frames.jl <outdir> [--seed 42] [--n 40] [--parcels 6] [--mcs 100] [--arms full,null]")
    outdir  = args[1]
    ispath(outdir) && error("destination already exists: $outdir")
    seed    = parse(Int, getopt(args, "--seed", "42"))
    N       = parse(Int, getopt(args, "--n", "40"))
    parcels = parse(Int, getopt(args, "--parcels", "6"))
    n_mcs   = parse(Int, getopt(args, "--mcs", "100"))
    arms    = String.(split(getopt(args, "--arms", "full,null"), ","))
    n_mcs > 0 || error("--mcs must be positive, got $n_mcs")
    allunique(arms) || error("duplicate arm in --arms")
    foreach(zeroes, arms)

    SR = load_serial()
    mkpath(outdir)
    frame0 = Dict{String, String}()
    counts = Dict{String, Vector{Int}}()
    for arm in arms
        p = arm_params(SR, arm, N, parcels)
        # One run per arm, replayed frame by frame: the state after k steps is what a run to
        # `k` would return, because the stream is seeded once and consumed in order. The
        # snapshot trajectory is not enough (it keeps totals, no lattice), so the lattice is
        # written from the state itself, re-run to each k. At N = 40 a run to 100 is seconds,
        # so 100 prefixes are minutes, and no step of the model is re-implemented here.
        counts[arm] = Int[]
        for k in 0:n_mcs
            state, _ = redirect_stdout(devnull) do
                Base.invokelatest(SR.run_simulation, p, k; seed = seed)
            end
            path = joinpath(outdir, @sprintf("%s_mcs%03d.csv", arm, k))
            push!(counts[arm], write_frame(path, SR, state))
            if k == 0
                frame0[arm] = read(path, String)
            end
            k % 10 == 0 && @printf(stderr, "%s mcs %d  %d sites\n", arm, k, counts[arm][end])
        end
    end
    # THE CONTROL: the seeding frame is one frame, whatever the arm.
    f0 = collect(values(frame0))
    all(==(f0[1]), f0) || error("the MCS 0 frames differ across arms; the arms do not share a placement")

    open(joinpath(outdir, "manifest.txt"), "w") do io
        println(io, "# consortium frames")
        println(io, "seed=$seed")
        println(io, "N=$N")
        println(io, "parcels_per_species=$parcels")
        println(io, "n_mcs=$n_mcs")
        println(io, "arms=", join(arms, ","))
        println(io, "arm_definitions=full: shipped CPMParams; no_rad: beta_ion = 0; no_mel: alpha_M_species = 0; null: both")
        println(io, "git_commit=", git_describe())
        println(io, "julia_version=$(VERSION)")
        println(io, "machine=$(Sys.MACHINE)")
        println(io, "species_names=", join(SR.SPECIES_NAMES, ";"))
        println(io, "frame_columns=x,y,z,species,cell")
        for arm in arms
            println(io, "occupied_sites_$arm=", join(counts[arm], " "))
        end
        println(io, "seeding_frame_identical_across_arms=true")
    end
    println("wrote $(length(arms) * (n_mcs + 1)) frames and manifest.txt to $outdir")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
