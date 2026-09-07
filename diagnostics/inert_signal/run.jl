#!/usr/bin/env julia
using HDF5, JSON3, SHA, TOML, Dates, WriteVTK
include("InertSignal.jl")
using .InertSignal
include(joinpath(@__DIR__, "..", "..", "export_vti.jl"))

sha(path) = bytes2hex(open(sha256, path))
require(ok, message) = ok || throw(ArgumentError(message))

function read_config(path)
    c = TOML.parsefile(path)
    expected = Set(["schema_version", "basis", "molecule", "production", "diffusion", "decay",
                    "spacing", "max_dt", "threshold", "source_timing", "physical_conversion",
                    "acceptance_coupling"])
    require(Set(keys(c)) == expected, "unknown or missing configuration key")
    require(c["schema_version"] == 1 && c["basis"] == "declared", "unrecognised declaration")
    require(c["molecule"] == "unassigned generic signal", "molecule-specific rates need an evidence audit")
    require(c["acceptance_coupling"] === false, "this diagnostic has no acceptance coupling")
    require(c["source_timing"] == "occupancy at MCS k held over [k,k+1)", "unknown source timing")
    require(c["physical_conversion"] == "blocked: D-PITCH and D-TIMESERIES", "physical conversion remains blocked")
    require(c["spacing"] == 1.0, "production diagnostic uses unit lattice spacing")
    p = SignalParams(; (Symbol(k) => c[k] for k in
        ("diffusion", "decay", "production", "spacing", "max_dt", "threshold"))...)
    c, p
end

"Validate the parent's pinned bytes before reading labels; no restart or uptake fields are read."
function verify_parent(parent, expected_sha)
    manifest_path = joinpath(parent, "run_manifest.json")
    require(sha(manifest_path) == expected_sha, "parent manifest hash mismatch")
    m = JSON3.read(read(manifest_path, String), Dict{String,Any})
    c = m["configuration"]
    require(c["N"] == 40 && c["seed"] == 42 && c["initial_parcels"] == 42 &&
            c["parcels_per_species"] == 6, "wrong manuscript configuration")
    rows = m["snapshots"]
    require([r["mcs"] for r in rows] == collect(0:100), "missing, duplicate or unordered MCS")
    files = ["snap_mcs$(lpad(t, 6, '0')).h5" for t in 0:100]
    require(sort(readdir(joinpath(parent, "snapshots"))) == files, "unexpected snapshot inventory")
    for (rel, hash) in m["artifacts"]
        require(!isabspath(rel) && !(".." in splitpath(rel)), "unsafe artifact path")
        require(sha(joinpath(parent, rel)) == hash, "parent artifact hash mismatch: $rel")
    end
    mask_ref = nothing
    for (t, row) in enumerate(rows)
        rel = "snapshots/" * files[t]
        require(row["path"] == rel && row["sha256"] == m["artifacts"][rel], "snapshot identity mismatch")
        h5open(joinpath(parent, rel), "r") do f
            a = attributes(f)
            require(read(a["mcs"]) == t-1, "snapshot MCS mismatch")
            require(read(a["logical_axis_order"]) == "xyz" &&
                    read(a["dataset_axis_order_h5py"]) == "zyx", "unknown snapshot axis order")
            require(read(a["cell_id_background"]) == 0 && read(a["cell_id_wall"]) == -1,
                    "unknown label sentinels")
            ids, sp = read(f["lattice/cell_id"]), read(f["lattice/species_id"])
            mask = Bool.(read(f["lattice/interior_mask"]))
            require(size(ids) == size(sp) == size(mask) == (40,40,40), "wrong grid shape")
            require(all(ids[.!mask] .== -1) && all(ids[mask] .>= 0), "mask/label inconsistency")
            require(mask_ref === nothing || mask == mask_ref, "mask changed")
            mask_ref = mask
            live, species = read(f["cells/id"]), read(f["cells/species"])
            present = sort(unique(ids[ids .> 0]))
            require(sort(live) == present && length(live) == 42, "registry/lattice parcel mismatch")
            require([count(==(s), species) for s in 1:7] == fill(6,7), "wrong per-species inventory")
            for (id,s) in zip(live, species)
                require(all(sp[ids .== id] .== s), "parcel species mismatch")
            end
            require(all(sp[(ids .== 0)] .== 0), "empty species mismatch")
            # The parent's hash includes its own explicit encoding; byte hashes above
            # are authoritative here, and its label receipt is carried without redefinition.
            require(read(a["label_state_hash"]) == row["label_state_hash"], "label receipt mismatch")
        end
    end
    m
end

function run_diagnostic(parent, output, parent_sha, config_path)
    require(!ispath(output), "output already exists; choose a new run directory")
    c, p = read_config(config_path)
    m = verify_parent(parent, parent_sha)
    root = normpath(joinpath(@__DIR__, "..", ".."))
    tracked_status = readchomp(`git -C $root status --porcelain --untracked-files=no`)
    require(isempty(tracked_status), "commit tracked source changes before production")
    commit = readchomp(`git -C $root rev-parse HEAD`)
    mkdir(output)
    mkpath(joinpath(output, "fields")); mkpath(joinpath(output, "paraview")); mkpath(joinpath(output, "provenance"))
    cp(config_path, joinpath(output, "provenance", "config.toml"))
    for name in ("Project.toml", "Manifest.toml", "InertSignal.jl", "run.jl")
        cp(joinpath(@__DIR__, name), joinpath(output, "provenance", name))
    end
    A = zeros(40,40,40)
    previous_species = nothing
    metrics = Any[]
    paraview_collection(joinpath(output, "paraview", "signal_trajectory")) do pvd
        for row in m["snapshots"]
            t = row["mcs"]
            src = joinpath(parent, row["path"])
            sp, mask = h5open(src, "r") do f
                read(f["lattice/species_id"]), Bool.(read(f["lattice/interior_mask"]))
            end
            balance = t == 0 ? nothing : update_signal!(A, previous_species, mask, p)
            previous_species = sp
            e = endpoint(A, sp, mask, p)
            push!(metrics, (; mcs=t, e..., maximum=maximum(A), mass_balance=balance))
            name = "signal_mcs$(lpad(t,6,'0'))"
            h5open(joinpath(output, "fields", name * ".h5"), "w") do f
                f["fields/signal", deflate=6] = A
                f["lattice/interior_mask", deflate=6] = UInt8.(mask)
                a = attributes(f)
                for (key, value) in ("mcs"=>t, "parent_manifest_sha256"=>parent_sha,
                    "parent_snapshot_sha256"=>row["sha256"], "parent_label_state_hash"=>row["label_state_hash"],
                    "logical_axis_order"=>"xyz", "dataset_axis_order_h5py"=>"zyx",
                    "signal_units"=>"declared arbitrary signal units", "time_units"=>"MCS",
                    "signal_sha256"=>bytes2hex(sha256(reinterpret(UInt8, vec(A)))),
                    "source_timing"=>c["source_timing"], "threshold"=>p.threshold,
                    "diffusion"=>p.diffusion, "decay"=>p.decay, "production"=>collect(p.production),
                    "dt_sub"=>(t == 0 ? 0.0 : balance.dt), "n_sub"=>(t == 0 ? 0 : balance.nsteps),
                    "wall_bc"=>"absorbing A=0 at masked walls", "cube_bc"=>"homogeneous Neumann by omitted flux",
                    "parameter_basis"=>"declared numerical demo", "spacing"=>p.spacing,
                    "acceptance_coupling"=>0, "c_s_analysis_blocked"=>1)
                    a[key] = value
                end
            end
            vtk = export_vti(src, joinpath(output, "paraview", name); save=false)
            vtk["signal", VTKCellData()] = A
            vtk["occupied_above_threshold", VTKCellData()] = UInt8.(mask .& (sp .> 0) .& (A .>= p.threshold))
            vtk["signal_units", VTKFieldData()] = "declared arbitrary units; inert generic signal; no molecule assigned"
            vtk["signal_source_timing", VTKFieldData()] = c["source_timing"]
            vtk["signal_threshold", VTKFieldData()] = p.threshold
            vtk["parent_manifest_sha256", VTKFieldData()] = parent_sha
            vtk_save(vtk); pvd[Float64(t)] = vtk
        end
    end
    write(joinpath(output, "metrics.json"), JSON3.write(metrics))
    source_hashes = Dict(f => sha(joinpath(root,f)) for f in split(readchomp(`git -C $root ls-files`), '\n') if isfile(joinpath(root,f)))
    artifacts = Dict(relpath(joinpath(dir,f),output) => sha(joinpath(dir,f)) for
                     (dir,_,files) in walkdir(output) for f in files)
    result = Dict("schema_version"=>1, "run_id"=>basename(output), "parent_run_id"=>m["run_id"],
        "parent_manifest_sha256"=>parent_sha, "configuration"=>c, "source_commit"=>commit,
        "source_hashes"=>source_hashes, "julia_version"=>string(VERSION),
        "julia_executable_sha256"=>sha(joinpath(Sys.BINDIR, Base.julia_exename())),
        "created_utc"=>string(now(UTC)), "rng"=>"none; deterministic field driven by pinned seed-42 labels",
        "snapshot_mcs"=>collect(0:100), "artifacts"=>artifacts,
        "scope"=>"inert generic signal on a fixed CPM trajectory; no biological QS response or physical calibration",
        "hash_contract"=>"all output files except derived_manifest.json and its sha256 receipt")
    path = joinpath(output, "derived_manifest.json")
    write(path, JSON3.write(result)); write(joinpath(output, "derived_manifest.sha256"), sha(path)*"  derived_manifest.json\n")
    verify_parent(parent, parent_sha) # original scientific evidence remains byte-identical
    println(JSON3.write((; output, parent_sha, frames=length(metrics), final=metrics[end])))
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 4 || error("usage: run.jl PARENT_RUN NEW_OUTPUT EXPECTED_PARENT_MANIFEST_SHA256 CONFIG_TOML")
    run_diagnostic(ARGS...)
end
