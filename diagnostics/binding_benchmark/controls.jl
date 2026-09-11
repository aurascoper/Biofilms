#!/usr/bin/env julia
# The control set for the closed diffusion-and-binding benchmark, on the frozen geometry.
#
#   julia --project=diagnostics/binding_benchmark diagnostics/binding_benchmark/controls.jl \
#         <parent_evidence_dir> <out_dir> [--config benchmark.toml]
#
# Runnable, not narrated. Each control states the question it asks, the number it
# measured, and whether it fired -- a control reported in prose is a claim, and a control
# that cannot fail is not a control. Two of these (K8, K12) are expected to be RED and
# are recorded as FIRES when they are; a green there would mean the ledger had quietly
# absorbed the defect it exists to expose.
using Dates, Printf
include(joinpath(@__DIR__, "setup.jl"))

with_params(p::Params; kw...) = Params(
    get(kw, :D_c, p.D_c), get(kw, :lambda, p.lambda), get(kw, :k_on, p.k_on),
    get(kw, :k_off, p.k_off), get(kw, :B0, p.B0), get(kw, :dB, p.dB),
    get(kw, :K, p.K), get(kw, :n, p.n), get(kw, :q_ext, p.q_ext))

"Run and sample. Returns the final state, the ledger, and the sampled series."
function series(geo, p, B, dt, nsteps, every; c0, b0 = 0.0, kwargs...)
    st = make_state(geo, B; c0 = c0, b0 = b0)
    led = Ledger(st, geo)
    out = Vector{NamedTuple}()
    run_to(st, geo, p, dt, nsteps, led;
           observe = (s, l, x) -> begin
               (s == 0 || s == nsteps || s % every == 0) || return
               Id, Ib, IB = inventory(x, geo)
               push!(out, (step = s, t = l.t, dissolved = Id, bound = Ib, capacity = IB))
           end, kwargs...)
    st, led, out
end

order(a, b) = log2(abs(a) / abs(b))

function controls(parent::String, config::String)
    cfg, p = read_config(config)
    geo, snap_path, snap_sha = frozen_geometry(parent, cfg)
    A = neighbour_fraction(geo.occupied)
    B = capacity_field(geo, p; A = A)
    dt = Float64(cfg["dt"]); total = Float64(cfg["total_time"])
    nsteps = round(Int, total / dt); every = max(1, nsteps ÷ 10)
    c0 = Float64(cfg["c0"]); b0 = Float64(cfg["b0"])
    results = Vector{Dict{String, Any}}()
    add!(name, question, verdict, detail) =
        push!(results, Dict{String, Any}("control" => name, "question" => question,
                                         "verdict" => verdict, "detail" => detail))

    # --- the production run, reused by several controls -------------------------------
    st, led, ser = series(geo, p, B, dt, nsteps, every; c0 = c0, b0 = b0)
    Id, Ib, IB = inventory(st, geo)
    scale = led.initial_dissolved + led.initial_bound

    # K1 -- matched total capacity ------------------------------------------------------
    U, Bbar = matched_uniform_capacity(B, geo.occupied)
    stu, ledu, seru = series(geo, p, U, dt, nsteps, every; c0 = c0, b0 = b0)
    cap_mismatch = maximum(abs(a.capacity - b.capacity) / a.capacity
                           for (a, b) in zip(ser, seru))
    dbound = maximum(abs(a.bound - b.bound) for (a, b) in zip(ser, seru))
    ddiss = maximum(abs(a.dissolved - b.dissolved) for (a, b) in zip(ser, seru))
    # Binding is internal and both pools decay at the same rate, so the total inventory
    # obeys dI/dt = -lambda I in both runs and must be identical between them at every
    # time. Any difference in the bound pool is therefore exactly mirrored in the
    # dissolved pool. This is free to check and it fails loudly if either run is wrong.
    dtotal = maximum(abs((a.dissolved + a.bound) - (b.dissolved + b.bound))
                     for (a, b) in zip(ser, seru))
    mirrored = abs(dbound - ddiss) <= 1e-9 * max(dbound, 1.0)
    add!("K1 matched total capacity",
         "at equal integrated capacity at every comparison time, does the spatial " *
         "arrangement of capacity change the trajectory?",
         (cap_mismatch < 1e-12 && dtotal < 1e-9 * scale && mirrored) ? "MEASURED" : "INVALID",
         Dict("uniform_capacity" => Bbar,
              "capacity_integral_local" => ser[end].capacity,
              "capacity_integral_uniform" => seru[end].capacity,
              "max_relative_capacity_mismatch" => cap_mismatch,
              "comparison_times" => [r.t for r in ser],
              "max_abs_bound_difference" => dbound,
              "max_relative_bound_difference" => dbound / maximum(r.bound for r in ser),
              "max_abs_dissolved_difference" => ddiss,
              "max_abs_total_inventory_difference" => dtotal,
              "bound_and_dissolved_differences_mirror" => mirrored,
              "note" => cap_mismatch < 1e-12 ?
                        "capacity integrals agree at every comparison time, so the " *
                        "difference below is attributable to arrangement, not amount" :
                        "capacity integrals differ; the comparison is not matched"))

    # K2 -- the dB = 0 reduction --------------------------------------------------------
    p0 = with_params(p; dB = 0.0)
    B0f = capacity_field(geo, p0; A = A)
    U0, _ = matched_uniform_capacity(B0f, geo.occupied)
    st0, _, ser0 = series(geo, p0, B0f, dt, nsteps, every; c0 = c0, b0 = b0)
    stu0, _, seru0 = series(geo, p0, U0, dt, nsteps, every; c0 = c0, b0 = b0)
    identical = (B0f == U0) && (st0.c == stu0.c) && (st0.b == stu0.b)
    add!("K2 dB = 0 reduction",
         "with the modulation switched off, do the local and uniform runs coincide " *
         "exactly rather than approximately?",
         identical ? "PASS" : "FAIL",
         Dict("capacity_fields_identical" => B0f == U0,
              "dissolved_identical" => st0.c == stu0.c,
              "bound_identical" => st0.b == stu0.b,
              "final_bound" => ser0[end].bound))

    # K3 -- the analytical binding limit -------------------------------------------------
    # c is clamped, so this is a different scheme -- an open one -- and the material the
    # clamp returns is booked as chemostat input. Capacity is uniform so every occupied
    # site follows the same scalar solution.
    pa = with_params(p; dB = 0.0, D_c = 0.0)
    Ba = capacity_field(geo, pa; A = A)
    q_bind = pa.k_on * c0 + pa.k_off + pa.lambda
    b_eq = pa.k_on * c0 * pa.B0 / q_bind
    exact = b_eq + (0.0 - b_eq) * exp(-q_bind * total)
    errs = Float64[]; dts = Float64[]
    for k in 0:2
        d = dt / 2^k
        sa, la, _ = series(geo, pa, Ba, d, round(Int, total / d), 10^9;
                           c0 = c0, b0 = 0.0, clamp_dissolved = true)
        site = findfirst(geo.occupied)
        push!(dts, d); push!(errs, abs(sa.b[site] - exact))
    end
    orders = [order(errs[i], errs[i + 1]) for i in 1:length(errs)-1]
    add!("K3 analytical binding limit",
         "does the integrated bound pool approach b_eq + (b0 - b_eq) exp(-q_bind t) at " *
         "first order in dt?",
         all(o -> abs(o - 1) < 0.1, orders) ? "PASS" : "FAIL",
         Dict("q_bind" => q_bind, "b_eq" => b_eq, "exact_at_T" => exact,
              "dt" => dts, "abs_error" => errs, "observed_order" => orders,
              "scheme_note" => "clamp_dissolved = true is an open system; the " *
                               "production run never uses it"))

    # K4 -- pure diffusion ----------------------------------------------------------------
    pd = with_params(p; lambda = 0.0, k_on = 0.0, k_off = 0.0)
    std_ = make_state(geo, zeros(size(B)); c0 = 0.0, b0 = 0.0)
    @inbounds for ix in eachindex(std_.c)
        geo.occupied[ix] && (std_.c[ix] = c0)      # a hot spot on the frozen labels
    end
    ledd = Ledger(std_, geo)
    spread0 = maximum(std_.c[geo.interior]) - minimum(std_.c[geo.interior])
    mx = maximum(std_.c[geo.interior]); mn = minimum(std_.c[geo.interior])
    monotone = true
    run_to(std_, geo, pd, dt, nsteps, ledd;
           observe = (s, l, x) -> begin
               s == 0 && return
               m1 = maximum(x.c[geo.interior]); m2 = minimum(x.c[geo.interior])
               (m1 <= mx + 1e-12 && m2 >= mn - 1e-12) || (monotone = false)
               mx = m1; mn = m2
           end)
    Idd, Ibd, _ = inventory(std_, geo)
    mean_c = Idd / (count(geo.interior) * prod(geo.spacing))
    add!("K4 pure diffusion limit",
         "with every reaction and loss rate zero, is the dissolved inventory conserved " *
         "and does the discrete maximum principle hold?",
         (abs(Idd - ledd.initial_dissolved) < 1e-10 * ledd.initial_dissolved && monotone &&
          Ibd == 0) ? "PASS" : "FAIL",
         Dict("initial_inventory" => ledd.initial_dissolved, "final_inventory" => Idd,
              "relative_change" => (Idd - ledd.initial_dissolved) / ledd.initial_dissolved,
              "bound_inventory" => Ibd,
              "maximum_principle_held" => monotone,
              "spread_initial" => spread0,
              "spread_final" => maximum(std_.c[geo.interior]) - minimum(std_.c[geo.interior]),
              "interior_mean" => mean_c,
              "transport_residual" => ledd.transport_residual,
              "closure_residual" => closure_residual(std_, geo, ledd)))

    # K5 -- pure decay ----------------------------------------------------------------------
    pl = with_params(p; D_c = 0.0, k_on = 0.0, k_off = 0.0)
    stl, ledl, _ = series(geo, pl, zeros(size(B)), dt, nsteps, 10^9; c0 = c0, b0 = 0.0)
    predicted = c0 * (1 - pl.lambda * dt)^nsteps
    devs = maximum(abs(stl.c[ix] - predicted) for ix in eachindex(stl.c) if geo.interior[ix])
    add!("K5 pure decay limit",
         "with transport and binding off, does the field follow the exact discrete " *
         "sequence c0 (1 - lambda dt)^n, and is that the same thing as being accurate?",
         devs < 1e-13 * c0 ? "PASS" : "FAIL",
         Dict("predicted_discrete" => predicted,
              "continuum" => c0 * exp(-pl.lambda * total),
              "max_abs_deviation_from_discrete" => devs,
              "discretisation_error_vs_continuum" =>
                  abs(predicted - c0 * exp(-pl.lambda * total)),
              "note" => "dt <= 1/lambda is a positivity condition. At dt = 1/lambda the " *
                        "scheme returns exactly zero against a continuum value of 1/e."))

    # K6 -- positivity and the capacity bound ------------------------------------------------
    overshoot_final = maximum(st.b[ix] - st.B[ix] for ix in eachindex(st.b) if geo.interior[ix])
    add!("K6 positivity and 0 <= b <= B",
         "over the production run, does either pool go negative, and does the bound " *
         "pool ever exceed the local capacity?",
         (led.min_dissolved >= 0 && led.min_bound >= 0 && overshoot_final <= 0) ?
             "PASS" : "FAIL",
         Dict("min_dissolved" => led.min_dissolved, "min_bound" => led.min_bound,
              "max_overshoot_before_release" => led.max_overshoot,
              "max_b_minus_B_final" => overshoot_final))

    # K7 / K8 -- the conservative release, and the same run with the transfer omitted -------
    half = nsteps ÷ 2
    stress = function (release::Bool)
        s = make_state(geo, B; c0 = c0, b0 = b0)
        l = Ledger(s, geo)
        run_to(s, geo, p, dt, half, l)
        s.B .*= 0.25                       # a forced capacity decrease, mid-run
        run_to(s, geo, p, dt, nsteps - half, l; release_to_dissolved = release)
        s, l
    end
    s7, l7 = stress(true)
    r7 = closure_residual(s7, geo, l7)
    over7 = maximum(s7.b[ix] - s7.B[ix] for ix in eachindex(s7.b) if geo.interior[ix])
    add!("K7 conservative release under a forced capacity decrease",
         "when capacity is cut mid-run, is the overflow moved to the dissolved pool " *
         "rather than created or destroyed?",
         (l7.released > 0 && abs(r7) < 1e-9 * scale && over7 <= 0) ? "FIRES" : "DID-NOT-FIRE",
         Dict("capacity_factor" => 0.25, "step_of_decrease" => half,
              "released" => l7.released, "dropped" => l7.dropped,
              "closure_residual" => r7, "closure_residual_relative" => r7 / scale,
              "max_b_minus_B_final" => over7))

    s8, l8 = stress(false)
    r8 = closure_residual(s8, geo, l8)
    add!("K8 the same transfer, deliberately omitted",
         "with the release credited nowhere, does the ledger report the loss instead of " *
         "absorbing it?",
         (l8.dropped > 0 && abs(r8) > 1e-6 * scale &&
          abs(r8 + l8.dropped) < 1e-9 * scale) ? "FIRES" : "DID-NOT-FIRE",
         Dict("released" => l8.released, "dropped" => l8.dropped,
              "closure_residual" => r8, "closure_residual_relative" => r8 / scale,
              "residual_plus_dropped" => r8 + l8.dropped,
              "expected" => "closure_residual == -dropped",
              "note" => "this is the red case. A PASS here would mean the identity had " *
                        "been written to include the term that records the defect."))

    # K9 -- timestep refinement -----------------------------------------------------------
    finals = Float64[]; rdts = Float64[]
    for k in 0:2
        d = dt / 2^k
        sk, _, _ = series(geo, p, B, d, round(Int, total / d), 10^9; c0 = c0, b0 = b0)
        _, Ibk, _ = inventory(sk, geo)
        push!(rdts, d); push!(finals, Ibk)
    end
    ord9 = order(finals[1] - finals[2], finals[2] - finals[3])
    add!("K9 timestep refinement",
         "does the bound inventory at T converge, and at what observed order?",
         abs(ord9 - 1) < 0.15 ? "PASS" : "FAIL",
         Dict("dt" => rdts, "final_bound" => finals,
              "differences" => [finals[1] - finals[2], finals[2] - finals[3]],
              "observed_order" => ord9))

    # K10 -- the two denominators ----------------------------------------------------------
    f_inv, f_cap = bound_fractions(st, geo)
    add!("K10 both bound-fraction denominators",
         "reported apart, because one number cannot answer both questions",
         "MEASURED",
         Dict("bound_fraction_of_inventory" => f_inv,
              "of_inventory_denominator" => Id + Ib,
              "bound_fraction_of_capacity" => f_cap,
              "of_capacity_denominator" => IB,
              "ratio" => f_cap / f_inv))

    # K11 -- ledger closure -----------------------------------------------------------------
    r11 = closure_residual(st, geo, led)
    # The bound pool closes on its own terms as well: everything that entered it came
    # through the reaction, and everything that left went to decay or to release. This is
    # a second identity over a subset of the same accumulators, so an error that happened
    # to cancel in the grand total still shows up here.
    r11b = Ib - (led.initial_bound + led.bound_transfer - led.decayed_from_bound -
                 led.released - led.dropped)
    add!("K11 ledger closure from independent accumulators",
         "does the inventory read off the fields agree with the inventory the step " *
         "terms booked, with nothing obtained by differencing?",
         (abs(r11) < 1e-10 * scale && abs(r11b) < 1e-10 * scale) ? "PASS" : "FAIL",
         Dict("initial_dissolved" => led.initial_dissolved,
              "initial_bound" => led.initial_bound,
              "final_dissolved" => Id, "final_bound" => Ib,
              "decayed_from_dissolved" => led.decayed_from_dissolved,
              "decayed_from_bound" => led.decayed_from_bound,
              "bound_transfer" => led.bound_transfer,
              "released" => led.released, "dropped" => led.dropped,
              "external_input" => led.external_input,
              "chemostat_input" => led.chemostat_input,
              "transport_residual" => led.transport_residual,
              "closure_residual" => r11, "closure_residual_relative" => r11 / scale,
              "bound_pool_residual" => r11b,
              "bound_pool_identity" =>
                  "I_b = I_b(0) + bound_transfer - decayed_from_bound - released - dropped"))

    # K12 -- the timestep restriction refuses what it should ---------------------------------
    rate0 = stability_rate(make_state(geo, B; c0 = c0, b0 = b0), geo, p)
    too_big = 1.0001 / rate0
    refused = false
    try
        sx = make_state(geo, B; c0 = c0, b0 = b0)
        step!(sx, geo, p, too_big, Ledger(sx, geo))
    catch e
        refused = e isa ErrorException
    end
    add!("K12 the derived restriction refuses a step above the bound",
         "is the bound enforced on the frozen geometry, or only documented?",
         refused ? "FIRES" : "DID-NOT-FIRE",
         Dict("rate" => rate0, "bound" => 1 / rate0, "attempted_dt" => too_big,
              "declared_dt" => dt, "declared_margin" => dt * rate0, "refused" => refused))

    cfg, geo, snap_path, snap_sha, results
end

function main(parent::String, out::String, config::String)
    refuse_existing(out)
    hashes = code_hashes()
    cfg, geo, snap_path, snap_sha, results = controls(parent, config)
    fresh_destination(out)
    doc = Dict{String, Any}(
        "diagnostic" => "closed diffusion-and-binding benchmark, control set",
        "configuration_file" => basename(config),
        "code_sha256" => hashes,
        "parent" => Dict("snapshot" => cfg["parent_snapshot"], "snapshot_sha256" => snap_sha,
                         "manifest_sha256" => cfg["parent_manifest_sha256"],
                         "run_id" => geo.run_id, "mcs" => geo.mcs),
        "expected_red" => ["K8 the same transfer, deliberately omitted"],
        "controls" => results,
        "julia_version" => string(VERSION),
        "created_utc" => string(now(UTC)) * "Z")
    write_json(joinpath(out, "control_verification.json"), doc)
    for r in results
        @printf("%-12s %-58s %s\n", r["verdict"], first(r["control"], 58), "")
    end
    bad = [r["control"] for r in results
           if r["verdict"] in ("FAIL", "INVALID", "DID-NOT-FIRE")]
    println(isempty(bad) ? "all controls reported as expected" :
            "unexpected verdicts: " * join(bad, "; "))
    println("wrote ", joinpath(out, "control_verification.json"))
    isempty(bad) || exit(1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: controls.jl <parent_evidence_dir> <out_dir> [--config <toml>]")
    i = findfirst(==("--config"), ARGS)
    config = isnothing(i) ? joinpath(@__DIR__, "benchmark.toml") : ARGS[i + 1]
    main(ARGS[1], ARGS[2], config)
end
