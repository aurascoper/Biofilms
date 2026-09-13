# Configuration reading and parent verification, shared by run.jl and controls.jl.
#
# The chain of custody has exactly one pin: `parent_manifest_sha256`. The manifest is
# hashed against it, the snapshot is hashed against the entry the manifest carries for
# the snapshot, and the snapshot's own declarations are then checked against the config.
# Carrying a `git_sha` would say which code produced the bytes; it would not say that
# these are the bytes.
using SHA, TOML, JSON3
include(joinpath(@__DIR__, "BindingBenchmark.jl"))
using .BindingBenchmark

sha256_file(path) = bytes2hex(open(SHA.sha256, path))
require(ok, msg) = ok || throw(ArgumentError(msg))

const CONFIG_KEYS = Set([
    "schema_version", "basis", "claim", "time_unit", "length_unit", "concentration_unit",
    "physical_conversion", "parent_manifest_sha256", "parent_snapshot", "parent_run_id",
    "parent_mcs", "parent_label_state_hash", "parent_mask_sha256", "capacity_modulator",
    "spacing", "dt", "total_time", "record_every", "c0", "b0", "params"])

const PARAM_KEYS = Set(["D_c", "lambda", "k_on", "k_off", "B0", "dB", "K", "n", "q_ext"])

"""
    read_config(path) -> (cfg, Params)

Refuse an unrecognised key rather than ignore it: a configuration this diagnostic
silently drops is a configuration a reader will believe was honoured.
"""
function read_config(path)
    cfg = TOML.parsefile(path)
    require(Set(keys(cfg)) == CONFIG_KEYS,
            "unknown or missing configuration key: " *
            "$(sort(collect(symdiff(Set(keys(cfg)), CONFIG_KEYS))))")
    require(Set(keys(cfg["params"])) == PARAM_KEYS,
            "unknown or missing parameter: " *
            "$(sort(collect(symdiff(Set(keys(cfg["params"])), PARAM_KEYS))))")
    require(cfg["schema_version"] == 1 && cfg["basis"] == "declared", "unrecognised declaration")
    require(cfg["physical_conversion"] == "blocked: D-PITCH and D-TIMESERIES",
            "physical conversion remains blocked for this diagnostic")
    require(cfg["time_unit"] == "dtu (declared diagnostic time unit)",
            "this benchmark runs on its own declared clock, not an imported one")
    require(cfg["capacity_modulator"] == "occupied face-neighbour fraction, 6-connected",
            "the capacity modulator is derived from the frozen labels by one declared rule")
    require(cfg["params"]["q_ext"] == 0.0, "a closed benchmark has no external source")
    require(cfg["spacing"] == [1.0, 1.0, 1.0],
            "lattice pitch is unit; a declared physical pitch is blocked")
    require(isfinite(cfg["dt"]) && cfg["dt"] > 0 && isfinite(cfg["total_time"]) && cfg["total_time"] > 0,
            "dt and total_time must be finite and positive")
    require(cfg["record_every"] isa Integer && cfg["record_every"] >= 1,
            "record_every must be a positive integer")
    # Domains. `--config` reaches this parser, so a sign or a NaN the shipped file never
    # carries can still arrive: a negative lambda is growth, a negative k_on is a sink at
    # zero, and either writes a receipt for a scheme this diagnostic does not describe.
    # TOML integers are widened here too, since Params is Float64-typed and `D_c = 0`
    # would otherwise be a MethodError rather than a configuration.
    p = Dict(k => Float64(v) for (k, v) in cfg["params"])
    for k in ("D_c", "lambda", "k_on", "k_off", "B0", "dB")
        require(isfinite(p[k]) && p[k] >= 0, "$k must be finite and non-negative, got $(p[k])")
    end
    require(isfinite(p["K"]) && p["K"] > 0 && isfinite(p["n"]) && p["n"] > 0,
            "Hill K and n must be finite and positive")
    require(p["B0"] + p["dB"] > 0, "capacity must be positive somewhere on the occupied set")
    c0, b0 = Float64(cfg["c0"]), Float64(cfg["b0"])
    require(isfinite(c0) && c0 >= 0 && isfinite(b0) && b0 >= 0,
            "c0 and b0 must be finite and non-negative")
    require(c0 + b0 > 0, "an empty initial inventory makes every relative residual 0/0")
    cfg, Params(p["D_c"], p["lambda"], p["k_on"], p["k_off"],
                p["B0"], p["dB"], p["K"], p["n"], p["q_ext"])
end

"""
    verify_parent(parent, cfg) -> (snapshot_path, snapshot_sha256)

Hash the manifest against the pin, then the snapshot against the manifest's own entry.
Every artifact path in the manifest is checked for traversal before it is joined.
"""
function verify_parent(parent::AbstractString, cfg)
    manifest_path = joinpath(parent, "run_manifest.json")
    isfile(manifest_path) || throw(ArgumentError("no run_manifest.json under $parent"))
    got = sha256_file(manifest_path)
    require(got == cfg["parent_manifest_sha256"],
            "parent manifest hash mismatch: $got != $(cfg["parent_manifest_sha256"])")
    m = JSON3.read(read(manifest_path, String), Dict{String, Any})
    require(m["run_id"] == cfg["parent_run_id"], "parent run_id mismatch")
    rel = cfg["parent_snapshot"]
    require(!isabspath(rel) && !(".." in splitpath(rel)), "unsafe snapshot path")
    expected = get(m["artifacts"], rel, nothing)
    require(!isnothing(expected), "the manifest does not carry $rel")
    path = joinpath(parent, rel)
    got_snap = sha256_file(path)
    require(got_snap == expected, "snapshot hash mismatch: $got_snap != $expected")
    path, got_snap
end

"""
    frozen_geometry(parent, cfg) -> (geo, snapshot_path, snapshot_sha256)

Verify, load, and cross-check. The snapshot's own `mcs`, `run_id`, `label_state_hash` and
`mask_sha256` must agree with the configuration: the byte hash proves the file was not
edited, and these prove it is the frame the configuration meant.
"""
function frozen_geometry(parent::AbstractString, cfg)
    path, snap_sha = verify_parent(parent, cfg)
    geo = load_geometry(path; spacing = Tuple(Float64.(cfg["spacing"])))
    require(geo.mcs == cfg["parent_mcs"], "snapshot is MCS $(geo.mcs), not $(cfg["parent_mcs"])")
    require(geo.run_id == cfg["parent_run_id"], "snapshot run_id mismatch")
    require(geo.label_state_hash == cfg["parent_label_state_hash"], "label state hash mismatch")
    require(geo.mask_sha256 == cfg["parent_mask_sha256"], "mask hash mismatch")
    require(count(geo.occupied) > 0, "the frozen frame carries no occupied sites")
    geo, path, snap_sha
end

"""
Steps in `total` at `dt`. A quotient that is not an integer is refused rather than
rounded: rounding integrates a different horizon while every control still compares
against the declared one. Shared by run.jl and controls.jl so they cannot disagree.
"""
function nsteps_for(total::Float64, dt::Float64)
    n = round(Int, total / dt)
    require(abs(n * dt - total) <= 1e-9 * max(1.0, total),
            "total_time / dt = $(total / dt) is not an integer number of steps")
    n
end

"Hashes of the code and configuration that produced a receipt, recorded beside it."
function code_hashes()
    Dict(f => sha256_file(joinpath(@__DIR__, f)) for f in
         ("BindingBenchmark.jl", "setup.jl", "run.jl", "controls.jl",
          "mutation_controls.jl", "test_numerics.jl", "benchmark.toml"))
end

"""
Refuse a destination that already exists. Checked twice on purpose: once before the run,
so a typo costs nothing, and once at write time, because a long run can finish into a
directory that appeared while it was running. The directory is created only at the
second call, so a run that fails leaves nothing behind to block the retry.

The second check is the creation itself. `mkpath` accepts a directory that already
exists, so a check followed by `mkpath` has a gap in which another process can create
the destination and have this run write into it; `mkdir` refuses with EEXIST in that
case, and nothing else can observe the gap. `race` is the test hook that reproduces it.
"""
refuse_existing(dir::AbstractString) =
    ispath(dir) && throw(ArgumentError("destination already exists: $dir"))

function fresh_destination(dir::AbstractString; race = () -> nothing)
    refuse_existing(dir)
    mkpath(dirname(abspath(dir)))
    race()
    mkdir(dir)
    dir
end

"""
The path as the filesystem would resolve it, for a path that need not exist yet: the
longest existing prefix is resolved through symlinks and the rest is appended lexically.
"""
function canonical(path::AbstractString)
    p = normpath(abspath(path)); rest = String[]
    while !ispath(p)
        pushfirst!(rest, basename(p)); p = dirname(p)
    end
    joinpath(realpath(p), rest...)
end

"""
Refuse an output at or below the parent bundle. `refuse_existing` cannot see this case:
a destination that does not exist yet passes it and is then created inside the evidence
the run promised not to write to.
"""
function refuse_inside(parent::AbstractString, out::AbstractString)
    p = canonical(parent); o = canonical(out)
    (o == p || startswith(o, p * "/")) &&
        throw(ArgumentError("output $out lies inside the parent bundle $parent"))
    nothing
end

write_json(path, obj) = open(path, "w") do io
    JSON3.pretty(io, JSON3.write(obj))
    println(io)
end
