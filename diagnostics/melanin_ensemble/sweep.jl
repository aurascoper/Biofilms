#!/usr/bin/env julia
# Seed sweep for the melanin observable, at a declared configuration.
#
#   julia diagnostics/melanin_ensemble/sweep.jl <out.csv> [--seeds 42:57]
#        [--n 20] [--parcels 2] [--mcs 400] [--at 100]
#
# This exists because the result it produces was previously run out of band. A post
# reporting "seeds 42 through 57" cited a commit whose tree contains no seed loop and no
# `--seed` flag: `biofilms_potts.jl` has exactly one `ARGS` reference and it is
# `--no-radiolysis`. A number nobody else can regenerate is a claim, not a measurement,
# and this repository already says so about controls.
#
# Nothing here re-implements the model. It loads `biofilms_potts.jl` through the same
# split-marker sandbox `validate_serial.jl` uses, so the sweep runs the shipped code and
# not a copy of it, and it emits the same CSV vocabulary: one row per (seed, mcs,
# species).
using Printf

"""
Load the serial model, minus its figure-export half.

The `#  13. Figure export` marker is load-bearing: seven files in this repository split
the monolith on that exact two-space string and evaluate only the half above it, in a
module whose imports are hardcoded to these four stdlibs. Loading this way means no
CairoMakie, which is why this sweep writes a CSV and leaves plotting to whoever wants it.
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

function main(args)
    isempty(args) && error("usage: sweep.jl <out.csv> [--seeds 42:57] [--n 20] [--parcels 2] [--mcs 400] [--at 100]")
    out = args[1]
    ispath(out) && error("destination already exists: $out")
    seeds   = parse_seeds(getopt(args, "--seeds", "42:57"))
    N       = parse(Int, getopt(args, "--n", "20"))
    parcels = parse(Int, getopt(args, "--parcels", "2"))
    n_mcs   = parse(Int, getopt(args, "--mcs", "400"))
    at      = parse(Int, getopt(args, "--at", "100"))

    SR = load_serial()
    # Snapshot cadence must land exactly on `at`, otherwise the row the analysis wants
    # does not exist. Asserted rather than assumed.
    interval = gcd(at, n_mcs)
    at % interval == 0 || error("snapshot interval $interval does not reach MCS $at")
    p = Base.invokelatest(SR.CPMParams; N = N, n_cells_per_species = parcels,
                          snapshot_interval = interval)
    alpha = Base.invokelatest(getfield, p, :α_M_species)
    names = SR.SPECIES_NAMES

    open(out, "w") do io
        println(io, "# melanin ensemble sweep")
        println(io, "# N=$N parcels_per_species=$parcels n_mcs=$n_mcs snapshot_interval=$interval")
        println(io, "# seeds=", join(seeds, " "))
        println(io, "# observable=mean melanin over OCCUPIED SITES (volume-weighted; ",
                    "biofilms_potts.jl take_snapshot). The mean of per-parcel means is a ",
                    "different quantity and is not computed here.")
        println(io, "# alpha_M is a declared input, not a measurement.")
        println(io, "seed,mcs,species,species_name,alpha_M,volume,n_cells,mean_melanin")
        for seed in seeds
            t0 = time()
            _, traj = redirect_stdout(devnull) do
                Base.invokelatest(SR.run_simulation, p, n_mcs; seed = seed)
            end
            any(s -> s.mcs == at, traj) || error("no snapshot at MCS $at for seed $seed")
            for snap in traj, sd in snap.species_data
                @printf(io, "%d,%d,%d,%s,%.4f,%d,%d,%.6f\n", seed, snap.mcs, sd.species,
                        replace(names[sd.species], "," => ""), alpha[sd.species],
                        sd.total_volume, sd.n_cells, sd.mean_melanin)
            end
            @printf(stderr, "seed %d  %.1f s\n", seed, time() - t0)
        end
    end
    println("wrote ", out)
end

# Parenthesised: `@__FILE__ && x` parses as `@__FILE__(&& x)`.
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
