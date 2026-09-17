#!/usr/bin/env julia
# Serial-reference CSV for validating biofilms_potts_jacc.jl.
#
# THIS CHECK HAS A DETECTION FLOOR, AND IT IS NOT REDUNDANT WITH THE ONE BESIDE
# IT. Perturbing the melanin coefficient in compute_delta_H_terms from 0.5 to
# 0.5000001 leaves every line below byte-identical; 0.5 to 50.0 fails it. A
# ~2e-7 shift in ΔH moves exp(-ΔH/T_cpm) by ~4e-8, and that has to straddle one
# of roughly 6.4M uniform draws before a single move flips and the trajectory
# diverges. So this file guards the TRAJECTORY, at the scale where a trajectory
# can be said to have changed -- it does not certify bit-exactness of the
# acceptance arithmetic.
#
# That is a third kind of check-that-cannot-fail: not always-green and not
# always-red, but CONDITIONALLY green with the condition unstated. Stated now.
# tests/delta_h_decomposition.jl is the exact half of the pair -- it asserts
# under `===` that the four terms sum to the scalar on every sampled move, with
# a control showing summation order is observable on 15.7% of them. Green here
# is necessary and not sufficient; neither file alone covers what both do.
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
    # basis_gate_ack: this harness reproduces a bit-for-bit CPM trajectory and
    # records NO radiodialysis quantity -- its CSV is CPM columns plus rd.m,
    # whose ODE has no X_total or X_red in it. The coupled loop reconstructs
    # RadiolysisParams with X_total = mean(compute_radial_biomass(...)), which
    # RADIODIALYSIS: BLOCKED gates, so stepping it needs this explicit
    # acknowledgement. It is not a claim that the basis is valid. Nothing here
    # may report c or s. tests/radiodialysis_basis_gate.jl holds this to it.
    rp = SerialRef.RadiolysisParams(Nr = 40, Ddot_R = 1.0, c_ext = 1.0,
                                    basis_gate_ack = true)
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
