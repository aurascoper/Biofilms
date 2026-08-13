#!/usr/bin/env julia
# Serial-reference CSV for validating biofilms_potts_jacc.jl.
# Runs the untouched serial script (minus its CairoMakie section) and
# prints the identical CSV lines:  CSV,seed,species,vol,ncells,mel,survived
#
#   julia validate_serial.jl [seeds...]   (default 42 43 44)

using Printf

function load_serial()
    src = read(joinpath(@__DIR__, "biofilms_potts.jl"), String)
    src = split(src, "#  13. Figure export")[1]
    SerialRef = Module(:SerialRef)
    Base.eval(SerialRef, :(using LinearAlgebra, Statistics, Random, Printf))
    Base.include_string(SerialRef, src, "biofilms_potts.jl")
    return SerialRef
end

function run(SerialRef::Module, seeds::Vector{Int})
    params = SerialRef.CPMParams(N = 40, n_cells_per_species = 6,
        snapshot_interval = 20)
    rp = SerialRef.RadiolysisParams(Nr = 40, Ddot_R = 1.0, c_ext = 1.0)
    for seed in seeds
        t0 = time()
        state, rd, _, _ = SerialRef.run_simulation_coupled(params, rp, 100; seed)
        snap = SerialRef.take_snapshot(state, 100)
        @printf("# seed %d done in %.1f s\n", seed, time() - t0)
        for sd in snap.species_data
            @printf("CSV,%d,%d,%d,%d,%.5f,%d\n", seed, sd.species,
                sd.total_volume, sd.n_cells, sd.mean_melanin,
                sd.n_cells > 0 ? 1 : 0)
        end
        @printf("CSVTOT,%d,%d,%.5f\n", seed, length(state.cells), rd.m)
    end
end

seeds = isempty(ARGS) ? [42, 43, 44] : parse.(Int, ARGS)
Base.invokelatest(run, load_serial(), seeds)
