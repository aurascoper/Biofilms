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
    require(cfg["dt"] > 0 && cfg["total_time"] > 0, "dt and total_time must be positive")
    require(cfg["record_every"] isa Integer && cfg["record_every"] >= 1,
            "record_every must be a positive integer")
    p = cfg["params"]
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
"""
refuse_existing(dir::AbstractString) =
    ispath(dir) && throw(ArgumentError("destination already exists: $dir"))

function fresh_destination(dir::AbstractString)
    refuse_existing(dir)
    mkpath(dir)
    dir
end

write_json(path, obj) = open(path, "w") do io
    JSON3.pretty(io, JSON3.write(obj))
    println(io)
end
