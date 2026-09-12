module LatticeEvidence

using HDF5, SHA, Serialization, TOML, Random, Statistics, Printf, Dates, JSON3
include("export_checkpoint.jl")
include("analysis/label_dynamics.jl")
const ROOT = @__DIR__
const N_MCS = 100
const CADENCES = (1, 2, 5, 10)
const GATE_NOTE = "RADIODIALYSIS: BLOCKED; c/s are retained only in complete restart checkpoints, never analysed or displayed by this diagnostic."

"""The manuscript coupled configuration; intentionally separate from the demo CLI.

basis_gate_ack is justified ONLY for labelled CPM diagnostics. The uptake
perturbation control exercises this factory and the production export. c/s
remain blocked, including at MCS 0; restarting requires keeping that contract.
"""
function manuscript_trajectory(SR)
    params = SR.CPMParams(N = 40, n_cells_per_species = 6, snapshot_interval = 20)
    rp = SR.RadiolysisParams(Nr = 40, Ddot_R = 1.0, c_ext = 1.0,
                            basis_gate_ack = true)
    return SR.init_coupled_simulation(params, rp; seed = 42)
end

filehash(p) = open(sha256, p) |> bytes2hex
arrayhash(a) = bytes2hex(sha256(reinterpret(UInt8, vec(Array(a)))))
params_dict(p) = Dict(string(k) => json_value(getfield(p, k)) for k in fieldnames(typeof(p)))
json_value(x::AbstractFloat) = isfinite(x) ? x : string(x)
json_value(x::AbstractArray) = map(json_value, x)
json_value(x) = x
read_json(p) = JSON3.read(read(p, String), Dict{String,Any})
function write_json(p, value)
    open(p, "w") do io
        JSON3.pretty(io, value)
        write(io, '\n')
    end
end

function label_arrays(SR, sim)
    species, lineage, generation = _label_arrays(SR, sim.state)
    return (; cell_id = Array{Int32}(sim.state.lattice), species, lineage, generation)
end
labelhash(a) = label_state_hash(a.cell_id, a.species, a.lineage, a.generation)

function inventory(ids, species, registry_ids, registry_species, volumes, mask)
    size(ids) == size(species) == size(mask) || error("grid mismatch")
    all(v -> v in (0, 1), mask) || error("invalid interior mask")
    inside = mask .== 1
    any(inside) || error("empty interior mask")
    all(ids[.!inside] .== -1) || error("wall sentinel mismatch")
    all(ids[inside] .>= 0) || error("negative interior label")
    all(species[ids .<= 0] .== 0) || error("nonzero background species")
    length(registry_ids) == length(unique(registry_ids)) || error("duplicate registry ID")
    length(registry_ids) == length(registry_species) == length(volumes) || error("registry shape mismatch")
    live = sort(Int.(registry_ids))
    lattice_ids = sort(Int.(unique(ids[ids .> 0])))
    live == lattice_ids || error("registry/lattice ID mismatch")
    all(s -> s in 1:7, registry_species) || error("unknown species")
    for (id, sp, vol) in zip(registry_ids, registry_species, volumes)
        sites = ids .== id
        count(sites) == vol > 0 || error("registry volume mismatch for ID $id")
        all(species[sites] .== sp) || error("parcel/species mismatch for ID $id")
    end
    counts = [count(==(s), registry_species) for s in 1:7]
    return Dict("live_registry_count" => length(live),
                "distinct_lattice_id_count" => length(lattice_ids),
                "per_species_parcel_counts" => counts,
                "per_species_site_counts" => [count(==(s), species[inside]) for s in 0:7],
                "interior_sites" => count(inside), "empty_interior_sites" => count(==(0), ids[inside]))
end

function inventory(SR, sim)
    a = label_arrays(SR, sim)
    cs = collect(values(sim.state.cells))
    return inventory(a.cell_id, a.species, collect(keys(sim.state.cells)),
                     [c.species for c in cs], [c.volume for c in cs], UInt8.(sim.state.interior))
end
function require_initial_inventory(info)
    info["live_registry_count"] == info["distinct_lattice_id_count"] == 42 || error("initial parcel count must be 42")
    info["per_species_parcel_counts"] == fill(6, 7) || error("initial per-species counts must all be six")
    return nothing
end

function write_snapshot(SR, sim, path, run_id; config = nothing)
    export_transport_snapshot(SR, sim, path; config_toml_path = config)
    info = inventory(SR, sim)
    a = label_arrays(SR, sim)
    h5open(path, "r+") do f
        attrs = attributes(f)
        attrs["run_id"] = run_id
        attrs["mask_sha256"] = arrayhash(UInt8.(sim.state.interior))
        attrs["grid_shape_xyz"] = collect(size(sim.state.lattice))
        attrs["live_registry_count"] = info["live_registry_count"]
        attrs["distinct_lattice_id_count"] = info["distinct_lattice_id_count"]
        attrs["per_species_parcel_counts"] = info["per_species_parcel_counts"]
        attrs["c_s_analysis_blocked"] = 1
        attrs["basis_gate_note"] = GATE_NOTE
    end
    merge!(info, Dict("mcs" => sim.mcs, "label_state_hash" => labelhash(a),
                     "sha256" => filehash(path)))
    return info
end

const INT_COLUMNS = (:mcs, :proposal_index, :donor_site, :recipient_site,
                     :donor_id, :recipient_id, :donor_species, :recipient_species)
const FLOAT_COLUMNS = (:adh, :vol, :rad, :mel, :delta_h, :draw)
mutable struct AcceptedCopies
    ints::Dict{Symbol,Vector{Int64}}
    floats::Dict{Symbol,Vector{Float64}}
end
AcceptedCopies() = AcceptedCopies(Dict(k => Int64[] for k in INT_COLUMNS),
                                  Dict(k => Float64[] for k in FLOAT_COLUMNS))
function (log::AcceptedCopies)(e)
    for k in INT_COLUMNS
        push!(log.ints[k], getproperty(e, k))
    end
    for k in FLOAT_COLUMNS
        push!(log.floats[k], getproperty(e, k))
    end
    nothing
end
function write_events(path, log, run_id)
    h5open(path, "w") do f
        a = attributes(f)
        a["schema_version"] = 1
        a["run_id"] = run_id
        a["logical_axis_order"] = "xyz"
        a["linear_index_base"] = 1
        a["linear_index_order"] = "Julia column-major; x fastest"
        a["proposal_index_semantics"] = "1-based index of all N^3 attempts in the MCS, including skips"
        a["draw_semantics"] = "NaN iff delta_h <= 0; no draw was made"
        a["observation_point"] = "after acceptance, before copy mutation"
        a["event_count"] = length(log.ints[:mcs])
        a["energy_terms"] = "adhesion, volume, radiation, melanin; acceptance units"
        for columns in (log.ints, log.floats), (k, v) in columns
            f["accepted/$(k)"] = v
        end
    end
end

function source_provenance(root = ROOT)
    tracked = filter(!isempty, split(read(`git -C $root ls-files -z`, String), '\0'))
    added = filter(!isempty, split(read(`git -C $root ls-files --others --exclude-standard -z`, String), '\0'))
    candidates = sort(unique(vcat(tracked, added)))
    hashes = Dict(p => filehash(joinpath(root, p)) for p in candidates if isfile(joinpath(root, p)))
    return Dict("commit" => readchomp(`git -C $root rev-parse HEAD`),
                "branch" => readchomp(`git -C $root branch --show-current`),
                "working_tree_status" => read(`git -C $root status --porcelain`, String),
                "tracked_source_hashes" => Dict(p => hashes[p] for p in tracked if haskey(hashes, p)),
                "additional_source_hashes" => Dict(p => hashes[p] for p in added if haskey(hashes, p)),
                "julia_version" => string(VERSION),
                "julia_threads" => Threads.nthreads(), "machine" => Sys.MACHINE,
                "julia_executable_sha256" => filehash(joinpath(Sys.BINDIR, Base.julia_exename())),
                "project_sha256" => filehash(joinpath(root, "Project.toml")),
                "manifest_sha256" => isfile(joinpath(root, "Manifest.toml")) ? filehash(joinpath(root, "Manifest.toml")) : nothing)
end

function refresh_artifacts!(dir, manifest)
    # A manifest cannot contain its own byte hash. Its separate receipt covers
    # it; all other run files are enumerated here, including restart files.
    files = Dict{String,String}()
    for (root, _, names) in walkdir(dir), name in names
        path = joinpath(root, name)
        rel = relpath(path, dir)
        rel in ("run_manifest.json", "run_manifest.sha256") && continue
        files[rel] = filehash(path)
    end
    manifest["artifacts"] = files
    manifest["hash_contract"] = "SHA-256 of every run file except this manifest and its separate run_manifest.sha256 receipt"
    write_json(joinpath(dir, "run_manifest.json"), manifest)
    write(joinpath(dir, "run_manifest.sha256"), filehash(joinpath(dir, "run_manifest.json")) * "\n")
    return manifest
end

function produce_run(SR, run_id; root = joinpath(ROOT, "out", "lattice_evidence"), config = nothing)
    occursin(r"^[A-Za-z0-9][A-Za-z0-9_-]*$", run_id) || error("run-id must be a safe directory name")
    config === nothing || isfile(config) || error("config file is absent")
    dir = joinpath(root, run_id)
    ispath(dir) && error("refusing existing run directory: $dir")
    sim = manuscript_trajectory(SR)
    initial = inventory(SR, sim)
    require_initial_inventory(initial)
    mkpath(root)
    mkdir(dir)                         # atomic refusal if another process wins
    mkdir(joinpath(dir, "snapshots"))
    mkdir(joinpath(dir, "restarts"))
    mkdir(joinpath(dir, "provenance"))
    for name in ("Project.toml", "Manifest.toml")
        src = joinpath(ROOT, name)
        isfile(src) && cp(src, joinpath(dir, "provenance", name))
    end
    config === nothing || cp(config, joinpath(dir, "provenance", "embedded_config.toml"))
    manifest = Dict{String,Any}("schema_version" => 1, "run_id" => run_id,
        "status" => "incomplete", "created_utc" => string(now(UTC)),
        "configuration" => Dict("N" => 40, "n_species" => 7, "initial_parcels" => 42,
            "parcels_per_species" => 6, "seed" => 42, "rng_type" => string(typeof(sim.rng)),
            "mcs_start" => 0, "mcs_end" => N_MCS, "cadence_mcs" => 1,
            "cpm_parameters" => params_dict(sim.state.params), "rd_parameters" => params_dict(sim.rd.params)),
        "provenance" => source_provenance(), "gate" => Dict("c_s_analysis_blocked" => true,
            "basis_gate_ack" => true, "note" => GATE_NOTE),
        "snapshots" => Any[], "initial_inventory" => initial,
        "figure1_hashes" => Dict(relpath(p, ROOT) => filehash(p) for p in
            filter(p -> startswith(basename(p), "fig1_radial_stratification"),
                   readdir(joinpath(ROOT, "preprint", "figures"); join = true))))
    write_json(joinpath(dir, "run_manifest.json"), manifest)
    export_restart_checkpoint(SR, sim, joinpath(dir, "restarts", "restart_mcs000000.h5"))
    # The underlying MCS-0 restart can truthfully declare a standalone basis.
    # This diagnostic's c/s prohibition is declared separately at every time.
    log = AcceptedCopies()
    for mcs in 0:N_MCS
        mcs > 0 && SR.advance_window!(sim, 1; on_accepted = log)
        path = joinpath(dir, "snapshots", @sprintf("snap_mcs%06d.h5", mcs))
        info = write_snapshot(SR, sim, path, run_id; config)
        info["path"] = relpath(path, dir)
        push!(manifest["snapshots"], info)
        mcs % 20 == 0 && @printf("labelled trajectory: MCS %d/%d\n", mcs, N_MCS)
    end
    export_restart_checkpoint(SR, sim, joinpath(dir, "restarts", "restart_mcs000100.h5"))
    for mcs in (0, N_MCS)
        h5open(joinpath(dir, "restarts", @sprintf("restart_mcs%06d.h5", mcs)), "r+") do f
            attributes(f)["run_id"] = run_id
            attributes(f)["c_s_analysis_blocked"] = 1
            attributes(f)["diagnostic_gate_note"] = GATE_NOTE
        end
    end
    write_events(joinpath(dir, "accepted_copies.h5"), log, run_id)
    manifest["accepted_copy_count"] = length(log.ints[:mcs])
    manifest["status"] = "exported_unverified"
    refresh_artifacts!(dir, manifest)
    return dir
end

function read_snapshot(path)
    h5open(path, "r") do f
        a = Dict(k => read(attributes(f)[k]) for k in keys(attributes(f)))
        arrays = (cell_id = read(f["lattice/cell_id"]), species = read(f["lattice/species_id"]),
                  lineage = read(f["lattice/lineage_id"]), generation = read(f["lattice/generation"]))
        mask = read(f["lattice/interior_mask"])
        info = inventory(arrays.cell_id, arrays.species, read(f["cells/id"]),
                         read(f["cells/species"]), read(f["cells/volume"]), mask)
        for key in ("live_registry_count", "distinct_lattice_id_count", "per_species_parcel_counts")
            info[key] == a[key] || error("snapshot $key mismatch: $path")
        end
        labelhash(arrays) == a["label_state_hash"] || error("altered snapshot label hash: $path")
        arrayhash(mask) == a["mask_sha256"] || error("altered mask hash: $path")
        collect(size(mask)) == a["grid_shape_xyz"] || error("grid metadata mismatch")
        a["logical_axis_order"] == "xyz" && a["dataset_axis_order_h5py"] == "zyx" || error("unknown axis order")
        a["cell_id_background"] == 0 && a["cell_id_wall"] == -1 || error("unknown sentinels")
        a["c_s_analysis_blocked"] == 1 || error("missing c/s analysis prohibition")
        return (; arrays..., mask, info, attrs = a)
    end
end

function verify_snapshots(dir; expected_mcs = collect(0:N_MCS))
    paths = filter(p -> endswith(p, ".h5"), readdir(joinpath(dir, "snapshots"); join = true))
    length(paths) == length(expected_mcs) || error("missing or extra snapshot MCS")
    frames = [read_snapshot(p) for p in paths]
    times = [Int(f.attrs["mcs"]) for f in frames]
    length(unique(times)) == length(times) || error("duplicate snapshot MCS")
    order = sortperm(times)
    times[order] == expected_mcs || error("missing or unexpected snapshot MCS")
    frames = frames[order]
    require_initial_inventory(first(frames).info)
    for f in frames
        f.mask == first(frames).mask || error("mask changed between snapshots")
        f.attrs["run_id"] == first(frames).attrs["run_id"] || error("mixed run IDs")
    end
    return frames
end

function verify_artifacts(dir, manifest)
    actual = Set{String}()
    for (root, _, names) in walkdir(dir), name in names
        rel = relpath(joinpath(root, name), dir)
        rel in ("run_manifest.json", "run_manifest.sha256") || push!(actual, rel)
    end
    actual == Set(keys(manifest["artifacts"])) || error("artifact inventory differs from manifest")
    for (rel, hash) in manifest["artifacts"]
        path = joinpath(dir, rel)
        isfile(path) || error("missing artifact: $rel")
        filehash(path) == hash || error("artifact hash mismatch: $rel")
    end
    filehash(joinpath(dir, "run_manifest.json")) == strip(read(joinpath(dir, "run_manifest.sha256"), String)) || error("manifest receipt mismatch")
end

function replay_events(dir, frames)
    rows, attrs = h5open(joinpath(dir, "accepted_copies.h5"), "r") do f
        Dict(k => read(f["accepted/$k"]) for k in keys(f["accepted"])),
        Dict(k => read(attributes(f)[k]) for k in keys(attributes(f)))
    end
    n = length(rows["mcs"])
    n == attrs["event_count"] > 0 || error("missing event record")
    all(length(v) == n for v in values(rows)) || error("event columns differ in length")
    attrs["run_id"] == first(frames).attrs["run_id"] || error("event run ID mismatch")
    attrs["linear_index_base"] == 1 && attrs["logical_axis_order"] == "xyz" || error("unknown event site indexing")
    ids = copy(frames[1].cell_id)
    species = copy(frames[1].species)
    lineage = copy(frames[1].lineage)
    generation = copy(frames[1].generation)
    mask = frames[1].mask .== 1
    last_order = (0, 0)
    cursor = 1
    counts_per_mcs = Int[]
    # Event-level return episodes are deliberately separate from sampled ones.
    ref_ids, ref_species = copy(ids), copy(species)
    event_returns = Dict("parcel" => zeros(Int32, size(ids)), "species" => zeros(Int32, size(ids)))
    for f in frames[2:end]
        mcs = Int(f.attrs["mcs"])
        before = cursor
        while cursor <= n && rows["mcs"][cursor] == mcs
            r(k) = rows[k][cursor]
            order = (Int(r("mcs")), Int(r("proposal_index")))
            order > last_order || error("events not strictly ordered")
            1 <= order[2] <= length(ids) || error("proposal index outside sweep")
            d, t = Int(r("donor_site")), Int(r("recipient_site"))
            1 <= d <= length(ids) && 1 <= t <= length(ids) || error("event site out of bounds")
            mask[d] && mask[t] || error("event outside interior")
            dc, tc = Tuple(CartesianIndices(ids)[d]), Tuple(CartesianIndices(ids)[t])
            maximum(abs.(collect(dc) .- collect(tc))) == 1 || error("non-neighbour copy")
            ids[d] == r("donor_id") && ids[t] == r("recipient_id") || error("pre-copy ID mismatch")
            species[d] == r("donor_species") && species[t] == r("recipient_species") || error("pre-copy species mismatch")
            ids[d] != ids[t] || error("same-ID copy in accepted log")
            dh = r("adh") + r("vol") + r("rad") + r("mel")
            isequal(dh, r("delta_h")) || error("energy sum mismatch")
            (dh <= 0 ? isnan(r("draw")) : 0 <= r("draw") < exp(-dh / 5.0)) || error("invalid Metropolis draw")
            event_returns["parcel"][t] += ids[t] != ref_ids[t] && ids[d] == ref_ids[t]
            event_returns["species"][t] += species[t] != ref_species[t] && species[d] == ref_species[t]
            ids[t], species[t], lineage[t], generation[t] = ids[d], species[d], lineage[d], generation[d]
            last_order = order
            cursor += 1
        end
        push!(counts_per_mcs, cursor - before)
        replay = (; cell_id = ids, species, lineage, generation)
        labelhash(replay) == f.attrs["label_state_hash"] || error("replay hash mismatch at MCS $mcs")
    end
    cursor == n + 1 || error("unconsumed or out-of-window events")
    return Dict("accepted_copy_count" => n, "counts_per_mcs" => counts_per_mcs,
                "replayed_mcs" => [Int(f.attrs["mcs"]) for f in frames],
                "event_return_episodes" => Dict(k => sum(v[mask]) for (k, v) in event_returns))
end

function analyse_run(dir, frames; verify_only = false)
    mkpath(joinpath(dir, "analysis"))
    mask = first(frames).mask .== 1
    species = [f.species for f in frames]
    parcels = [f.cell_id for f in frames]
    summary = Dict{String,Any}("run_id" => first(frames).attrs["run_id"],
        "reference_mcs" => 0, "last_mcs" => N_MCS,
        "semantics" => "fixed interior; empty label 0 included; sampled observations, not within-MCS events",
        "interior_sites" => count(mask), "cadences" => Dict{String,Any}())
    h5open(joinpath(dir, "analysis", "label_dynamics.h5"), verify_only ? "r" : "w") do f
        function output(path, value)
            if verify_only
                haskey(f, path) && isequal(read(f[path]), value) || error("analysis mismatch: $path")
            else
                f[path] = value
            end
        end
        if verify_only
            read(attributes(f)["run_id"]) == summary["run_id"] || error("analysis run ID mismatch")
        else
            attributes(f)["run_id"] = summary["run_id"]
            attributes(f)["logical_axis_order"] = "xyz"
            attributes(f)["map_outside_domain"] = "counts zero; exclude with interior_mask"
            attributes(f)["occupancy_outside_domain"] = "NaN"
        end
        output("interior_mask", UInt8.(mask))
        output("initial_species", species[1])
        output("initial_parcel_id", parcels[1])
        for stride in CADENCES
            row = Dict{String,Any}("cadence_mcs" => stride, "mcs" => collect(0:stride:N_MCS))
            output("cadence_$stride/species_occupancy", LabelDynamics.species_occupancy(species, mask; stride))
            for (kind, labels) in (("species", species), ("parcel", parcels))
                metrics = LabelDynamics.sampled_history(labels, mask, species[1]; stride)
                prefix = "cadence_$stride/$kind"
                output("$prefix/transitions", metrics.transitions)
                output("$prefix/returns", metrics.returns)
                output("$prefix/hidden_reversal_intervals", metrics.hidden)
                output("$prefix/persistent", UInt8.(metrics.persistent))
                output("$prefix/endpoint_identity", UInt8.(metrics.endpoint_identity))
                row[kind] = metrics.rows
            end
            summary["cadences"][string(stride)] = row
        end
    end
    summary_path = joinpath(dir, "analysis", "label_dynamics.json")
    if verify_only
        read_json(summary_path) == summary || error("analysis summary mismatch")
    else
        write_json(summary_path, summary)
    end
    return summary
end

# Recursive equality includes every parameter, field, registry/event entry and
# RNG bytes. HDF5 container bytes need not match for equal scientific state.
function complete_equal(a, b)
    typeof(a) == typeof(b) || return false
    if a isa AbstractRNG
        ioa, iob = IOBuffer(), IOBuffer()
        serialize(ioa, a); serialize(iob, b)
        return take!(ioa) == take!(iob)
    elseif a isa AbstractDict
        Set(keys(a)) == Set(keys(b)) || return false
        return all(complete_equal(a[k], b[k]) for k in keys(a))
    elseif a isa AbstractArray
        axes(a) == axes(b) || return false
        return all(complete_equal(x, y) for (x, y) in zip(a, b))
    elseif a isa Number || a isa Symbol || a isa AbstractString || a === nothing
        return isequal(a, b)
    end
    return all(complete_equal(getfield(a, k), getfield(b, k)) for k in fieldnames(typeof(a)))
end

function perturb_uptake!(SR, sim)
    p = sim.rd.params
    sim.rd.params = SR.RadiolysisParams(;
        (k => getfield(p, k) for k in fieldnames(typeof(p)) if k ∉ (:k_ads, :k_red))...,
        k_ads = 10p.k_ads, k_red = 10p.k_red)
    return sim
end

function verify_determinism(SR, dir, frames)
    initial = restore_restart_checkpoint(SR, joinpath(dir, "restarts", "restart_mcs000000.h5"))
    finish = restore_restart_checkpoint(SR, joinpath(dir, "restarts", "restart_mcs000100.h5"))
    complete_equal(initial, manuscript_trajectory(SR)) || error("initial restart differs from canonical factory")
    one_window = manuscript_trajectory(SR)
    SR.advance_window!(one_window, N_MCS)
    complete_equal(finish, one_window) || error("instrumentation/windowing changed complete final state")
    windows = manuscript_trajectory(SR)
    perturbed = perturb_uptake!(SR, manuscript_trajectory(SR))
    for mcs in 1:N_MCS
        SR.advance_window!(windows, 1)
        SR.advance_window!(perturbed, 1)
        target = frames[mcs + 1].attrs["label_state_hash"]
        labelhash(label_arrays(SR, windows)) == target || error("uninstrumented label mismatch at MCS $mcs")
        labelhash(label_arrays(SR, perturbed)) == target || error("labelled outputs depend on gated uptake at MCS $mcs")
    end
    complete_equal(windows, one_window) || error("one-MCS and 100-MCS complete states differ")
    # Confirm continuation, not merely current labels or a displayed RNG name.
    for sim in (finish, one_window, windows)
        SR.advance_window!(sim, 1)
    end
    complete_equal(finish, one_window) && complete_equal(windows, one_window) || error("RNG continuation differs")
    shared = collect(0:20:N_MCS)
    for mcs in shared
        baseline = manuscript_trajectory(SR)
        result = redirect_stdout(devnull) do
            SR.run_simulation_coupled(baseline.state.params, baseline.rd.params, mcs; seed = 42)
        end
        state = result[1]
        species, lineage, generation = _label_arrays(SR, state)
        hash = label_state_hash(state.lattice, species, lineage, generation)
        hash == frames[mcs + 1].attrs["label_state_hash"] || error("legacy coupled-run mismatch at MCS $mcs")
    end
    return Dict("instrumented_uninstrumented_complete_state" => true,
                "one_and_100_mcs_windows_complete_state" => true,
                "rng_continuation_mcs" => 101, "legacy_shared_mcs" => shared,
                "uptake_perturbation_factor" => 10, "uptake_control_matched_mcs" => collect(1:N_MCS))
end

function postflight(SR, dir; determinism = true)
    manifest = read_json(joinpath(dir, "run_manifest.json"))
    verify_artifacts(dir, manifest)
    frames = verify_snapshots(dir)
    manifest["run_id"] == first(frames).attrs["run_id"] || error("manifest run ID mismatch")
    length(manifest["snapshots"]) == length(frames) || error("manifest snapshot count mismatch")
    for (record, frame) in zip(manifest["snapshots"], frames)
        record["mcs"] == frame.attrs["mcs"] || error("manifest MCS mismatch")
        record["label_state_hash"] == frame.attrs["label_state_hash"] || error("manifest label hash mismatch")
        record["sha256"] == filehash(joinpath(dir, record["path"])) || error("manifest snapshot hash mismatch")
        for key in keys(frame.info)
            record[key] == frame.info[key] || error("manifest snapshot $key mismatch")
        end
    end
    replay = replay_events(dir, frames)
    analysis = analyse_run(dir, frames; verify_only = haskey(manifest["artifacts"], "analysis/label_dynamics.h5"))
    verification = Dict{String,Any}("snapshot_count" => length(frames), "replay" => replay,
        "determinism_executed" => determinism)
    determinism && (verification["determinism"] = verify_determinism(SR, dir, frames))
    for (path, hash) in manifest["figure1_hashes"]
        filehash(joinpath(ROOT, path)) == hash || error("Figure 1 changed")
    end
    write_json(joinpath(dir, "verification.json"), verification)
    manifest["status"] = determinism ? (manifest["status"] == "manuscript_built" ? "manuscript_built" : "trajectory_verified") : "replay_verified_only"
    manifest["verification"] = verification
    refresh_artifacts!(dir, manifest)
    return (; manifest, analysis)
end

function main(args)
    length(args) >= 2 || error("usage: lattice_evidence.jl <run RUN_ID [--config FILE] | verify RUN_DIR>")
    SR = load_serial()
    Base.invokelatest() do
        if args[1] == "run"
            length(args) in (2, 4) || error("unrecognised runner arguments")
            cfg = length(args) == 4 ? (args[3] == "--config" ? args[4] : error("only --config is accepted")) : nothing
            dir = produce_run(SR, args[2]; config = cfg)
            postflight(SR, dir)
            println("verified trajectory: ", dir)
        elseif args[1] == "verify" && length(args) == 2
            postflight(SR, abspath(args[2]))
            println("postflight passed: ", abspath(args[2]))
        else
            error("unknown operation or arguments")
        end
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    LatticeEvidence.main(ARGS)
end
