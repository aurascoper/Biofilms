#!/usr/bin/env julia
# The closed diffusion-and-binding benchmark, on frozen MCS-0 labels and mask.
#
#   julia --project=diagnostics/binding_benchmark diagnostics/binding_benchmark/run.jl \
#         <parent_evidence_dir> <out_dir> [--config benchmark.toml]
#
# Writes `benchmark_receipt.json` and `timeseries.csv` into a destination that must not
# already exist. Nothing in the parent bundle is written to, and no figure, CSV or
# manuscript artifact anywhere in the repository is touched.
using Dates, Printf
include(joinpath(@__DIR__, "setup.jl"))

function modulator_histogram(A, occupied)
    counts = zeros(Int, 7)
    @inbounds for ix in eachindex(A)
        occupied[ix] && (counts[round(Int, A[ix] * 6) + 1] += 1)
    end
    Dict("levels" => ["$(i)/6" for i in 0:6], "counts" => counts)
end

function main(parent::String, out::String, config::String)
    refuse_existing(out)
    hashes = code_hashes()
    cfg, p = read_config(config)
    geo, snap_path, snap_sha = frozen_geometry(parent, cfg)

    A = neighbour_fraction(geo.occupied)
    B = capacity_field(geo, p; A = A)
    st = make_state(geo, B; c0 = Float64(cfg["c0"]), b0 = Float64(cfg["b0"]))
    led = Ledger(st, geo)

    # The restriction is derived from this stencil and these rates, then compared with
    # the declared step. It is re-evaluated inside every step as well, because max(c) and
    # max(B - b) both move and a run can leave the admissible region after entering it.
    dt = Float64(cfg["dt"])
    rate0 = stability_rate(st, geo, p)
    bound = 1 / rate0
    require(dt <= bound, "declared dt = $dt exceeds the derived bound $bound")

    total = Float64(cfg["total_time"])
    nsteps = round(Int, total / dt)
    require(abs(nsteps * dt - total) <= 1e-9 * max(1.0, total),
            "total_time / dt = $(total / dt) is not an integer number of steps")
    every = Int(cfg["record_every"])

    rows = Vector{Dict{String, Any}}()
    function observe(step, l, s)
        (step == 0 || step == nsteps || step % every == 0) || return
        Id, Ib, IB = inventory(s, geo)
        f_inv, f_cap = bound_fractions(s, geo)
        push!(rows, Dict{String, Any}(
            "step" => step, "t" => l.t,
            "dissolved" => Id, "bound" => Ib, "capacity" => IB,
            "bound_fraction_of_inventory" => f_inv,
            "bound_fraction_of_capacity" => f_cap,
            "decayed_from_dissolved" => l.decayed_from_dissolved,
            "decayed_from_bound" => l.decayed_from_bound,
            "bound_transfer" => l.bound_transfer,
            "released" => l.released,
            "transport_residual" => l.transport_residual,
            "closure_residual" => closure_residual(s, geo, l),
            "min_dissolved" => l.min_dissolved,
            "min_bound" => l.min_bound,
            "max_overshoot" => l.max_overshoot,
            "stability_margin" => dt * stability_rate(s, geo, p)))
    end

    elapsed = @elapsed run_to(st, geo, p, dt, nsteps, led; observe = observe)
    require(abs(led.t - nsteps * dt) <= 1e-9 * max(1.0, total),
            "the integrator's own clock drifted from steps x dt")

    Id, Ib, IB = inventory(st, geo)
    resid = closure_residual(st, geo, led)
    scale = led.initial_dissolved + led.initial_bound

    fresh_destination(out)
    open(joinpath(out, "timeseries.csv"), "w") do io
        cols = ["step", "t", "dissolved", "bound", "capacity",
                "bound_fraction_of_inventory", "bound_fraction_of_capacity",
                "decayed_from_dissolved", "decayed_from_bound", "bound_transfer",
                "released", "transport_residual", "closure_residual",
                "min_dissolved", "min_bound", "max_overshoot", "stability_margin"]
        println(io, join(cols, ","))
        for r in rows
            println(io, join((r[c] isa Integer ? string(r[c]) : @sprintf("%.17g", r[c])
                              for c in cols), ","))
        end
    end

    receipt = Dict{String, Any}(
        "diagnostic" => "closed diffusion-and-binding benchmark, Workstream B Stage 1",
        "claim" => cfg["claim"],
        "configuration" => cfg,
        "code_sha256" => hashes,
        # Paths are recorded as the configuration gave them, not as absolute paths on
        # the machine that ran this: the hashes are what bind the receipt to the bytes,
        # and an absolute path under one home directory binds nothing in a clone.
        "parent" => Dict("snapshot" => cfg["parent_snapshot"],
                         "snapshot_sha256" => snap_sha,
                         "manifest_sha256" => cfg["parent_manifest_sha256"],
                         "run_id" => geo.run_id, "mcs" => geo.mcs,
                         "label_state_hash" => geo.label_state_hash,
                         "mask_sha256" => geo.mask_sha256),
        "geometry" => Dict("grid" => collect(size(geo.occupied)),
                           "spacing" => collect(geo.spacing),
                           "interior_sites" => count(geo.interior),
                           "occupied_sites" => count(geo.occupied),
                           "modulator_histogram" => modulator_histogram(A, geo.occupied),
                           "modulator_min" => minimum(A[geo.occupied]),
                           "modulator_max" => maximum(A[geo.occupied]),
                           "capacity_min" => minimum(B[geo.occupied]),
                           "capacity_max" => maximum(B[geo.occupied]),
                           "capacity_integral" => IB),
        "timestep" => Dict(
            "rule" => "dt * max(2 D sum_j h_j^-2 + k_on max(B-b) + lambda, " *
                      "k_on max(c) + k_off + lambda) <= 1",
            "derived_from" => "the implemented stencil and the unsplit right-hand side",
            "diffusion_rate" => diffusion_rate(p.D_c, geo.spacing),
            "initial_rate" => rate0, "bound" => bound, "declared_dt" => dt,
            "initial_margin" => dt * rate0,
            "max_margin_observed" => maximum(r["stability_margin"] for r in rows),
            "steps" => nsteps),
        "final" => Dict(
            "t" => led.t, "dissolved" => Id, "bound" => Ib, "capacity" => IB,
            "bound_fraction_of_inventory" => bound_fractions(st, geo)[1],
            "bound_fraction_of_capacity" => bound_fractions(st, geo)[2],
            "initial_dissolved" => led.initial_dissolved,
            "initial_bound" => led.initial_bound,
            "decayed_from_dissolved" => led.decayed_from_dissolved,
            "decayed_from_bound" => led.decayed_from_bound,
            "bound_transfer" => led.bound_transfer,
            "released" => led.released, "dropped" => led.dropped,
            "external_input" => led.external_input,
            "chemostat_input" => led.chemostat_input,
            "transport_residual" => led.transport_residual,
            "closure_residual" => resid,
            "closure_residual_relative" => resid / scale,
            "min_dissolved" => led.min_dissolved, "min_bound" => led.min_bound,
            "max_overshoot" => led.max_overshoot),
        "timeseries" => rows,
        "elapsed_seconds" => elapsed,
        "julia_version" => string(VERSION),
        "created_utc" => string(now(UTC)) * "Z")
    write_json(joinpath(out, "benchmark_receipt.json"), receipt)

    @printf("interior %d, occupied %d, capacity integral %.6f\n",
            count(geo.interior), count(geo.occupied), IB)
    @printf("dt %.6g, derived bound %.6g, margin %.4f over %d steps (%.2f s)\n",
            dt, bound, dt * rate0, nsteps, elapsed)
    @printf("dissolved %.6f -> %.6f, bound %.6f -> %.6f\n",
            led.initial_dissolved, Id, led.initial_bound, Ib)
    @printf("decayed %.6f (%.6f dissolved + %.6f bound), released %.6g, dropped %.6g\n",
            led.decayed_from_dissolved + led.decayed_from_bound,
            led.decayed_from_dissolved, led.decayed_from_bound, led.released, led.dropped)
    @printf("bound fraction: %.6f of inventory, %.6f of capacity\n",
            bound_fractions(st, geo)...)
    @printf("closure residual %.6g (%.3g relative), transport residual %.6g\n",
            resid, resid / scale, led.transport_residual)
    @printf("min dissolved %.6g, min bound %.6g, max overshoot %.6g\n",
            led.min_dissolved, led.min_bound, led.max_overshoot)
    println("wrote ", joinpath(out, "benchmark_receipt.json"))
    receipt
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: run.jl <parent_evidence_dir> <out_dir> [--config <toml>]")
    i = findfirst(==("--config"), ARGS)
    config = isnothing(i) ? joinpath(@__DIR__, "benchmark.toml") : ARGS[i + 1]
    main(ARGS[1], ARGS[2], config)
end
