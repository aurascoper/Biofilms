#!/usr/bin/env julia
# HDF5 exchange writer/reader for the CPM ↔ OpenMC coupling.
# Two distinct schemas (docs/exchange_schema.md):
#   transport snapshot  — portable labels/fields for the Python/OpenMC side
#   restart checkpoint  — complete Julia state, bit-exact resume (RNG included)
#
# CLI (demo runs; tests call the functions directly):
#   julia --project=. export_checkpoint.jl transport out.h5 [--config cfg.toml] [--seed 42] [--mcs 20]
#   julia --project=. export_checkpoint.jl restart  out.h5 [--seed 42] [--mcs 20]

using HDF5, SHA, Serialization, TOML, Random, Statistics, Printf

# Same loading trick as validate_serial.jl / tests: monolith minus figures.
function load_serial()
    src = read(joinpath(@__DIR__, "biofilms_potts.jl"), String)
    src = split(src, "#  13. Figure export")[1]
    M = Module(:SerialRef)
    Base.eval(M, :(using LinearAlgebra, Statistics, Random, Printf))
    Base.include_string(M, src, "biofilms_potts.jl")
    return M
end

_git_sha() = try
    readchomp(`git -C $(@__DIR__) rev-parse HEAD`)
catch
    "unknown"
end

function _common_attrs!(f, state)
    a = attributes(f)
    a["schema_version"] = 1
    a["logical_axis_order"] = "xyz"
    a["dataset_axis_order_h5py"] = "zyx"
    a["coordinate_index_base"] = 0
    a["cell_id_background"] = 0
    a["cell_id_wall"] = -1
    a["git_sha"] = _git_sha()
    a["julia_version"] = string(VERSION)
    a["mcs"] = state.current_mcs
    a["physical_time_s"] = state.physical_time_s
    return f
end

# Per-site label arrays derived from the lattice + live-cell registry.
function _label_arrays(SR, state)
    N = state.params.N
    species = zeros(Int32, N, N, N)
    lineage = zeros(Int32, N, N, N)
    generation = zeros(Int32, N, N, N)
    for z in 1:N, y in 1:N, x in 1:N
        σ = state.lattice[x, y, z]
        σ > 0 || continue
        c = get(state.cells, Int(σ), nothing)
        c === nothing && continue
        species[x, y, z] = Int32(c.species)
        lineage[x, y, z] = Int32(c.lineage_id)
        generation[x, y, z] = Int32(c.generation)
    end
    return species, lineage, generation
end

label_state_hash(cell_id, species, lineage, generation) =
    bytes2hex(sha256(vcat(reinterpret(UInt8, vec(cell_id)),
                          reinterpret(UInt8, vec(species)),
                          reinterpret(UInt8, vec(lineage)),
                          reinterpret(UInt8, vec(generation)))))

# ≥3 asymmetric probe points (0-based), each with a distinct coordinate
# pattern so any axis swap or flip changes at least one expected value.
function _orientation_probes(cell_id)
    N = size(cell_id, 1)
    pts0 = [(1, 2, 3), (N - 3, 1, 0), (0, N - 4, N - 2), (2, N - 2, 1)]
    probes = Matrix{Int64}(undef, length(pts0), 4)
    for (i, (x, y, z)) in enumerate(pts0)
        probes[i, :] .= (x, y, z, Int64(cell_id[x + 1, y + 1, z + 1]))
    end
    return probes
end

function _write_cells!(f, state)
    ids = sort!(collect(keys(state.cells)))
    cs = [state.cells[i] for i in ids]
    f["cells/id"] = Int64.(ids)
    f["cells/species"] = Int64[c.species for c in cs]
    f["cells/volume"] = Int64[c.volume for c in cs]
    f["cells/lineage_id"] = Int64[c.lineage_id for c in cs]
    f["cells/parent_id"] = Int64[c.parent_id for c in cs]
    f["cells/generation"] = Int64[c.generation for c in cs]
    f["cells/birth_mcs"] = Int64[c.birth_mcs for c in cs]
    f["cells/lifetime_dose_Gy"] = Float64[c.lifetime_dose_Gy for c in cs]
    com = Matrix{Float64}(undef, 3, length(cs))
    for (j, c) in enumerate(cs)
        com[:, j] .= c.com
    end
    f["cells/com"] = com
end

function _write_events!(f, state)
    dl = state.division_log
    for (name, get_) in (
            ("mcs", e -> e.mcs), ("physical_time_s", e -> e.physical_time_s),
            ("parent_id", e -> e.parent_id),
            ("daughter_a_id", e -> e.daughter_a_id),
            ("daughter_b_id", e -> e.daughter_b_id),
            ("lineage_id", e -> e.lineage_id),
            ("parent_generation", e -> e.parent_generation),
            ("daughter_generation", e -> e.daughter_generation),
            ("parent_volume", e -> e.parent_volume),
            ("daughter_a_volume", e -> e.daughter_a_volume),
            ("daughter_b_volume", e -> e.daughter_b_volume))
        v = [get_(e) for e in dl]
        f["events/division/$name"] = name == "physical_time_s" ?
            Float64.(v) : Int64.(v)
    end
    de = state.death_log
    for (name, get_) in (
            ("mcs", e -> e.mcs), ("physical_time_s", e -> e.physical_time_s),
            ("cell_id", e -> e.cell_id), ("species", e -> e.species),
            ("lineage_id", e -> e.lineage_id), ("generation", e -> e.generation),
            ("lifetime_dose_Gy", e -> e.lifetime_dose_Gy))
        v = [get_(e) for e in de]
        f["events/death/$name"] = name in ("physical_time_s", "lifetime_dose_Gy") ?
            Float64.(v) : Int64.(v)
    end
    ar = state.archive
    f["archive/id"] = Int64[a.id for a in ar]
    f["archive/species"] = Int64[a.species for a in ar]
    f["archive/lineage_id"] = Int64[a.lineage_id for a in ar]
    f["archive/parent_id"] = Int64[a.parent_id for a in ar]
    f["archive/generation"] = Int64[a.generation for a in ar]
    f["archive/birth_mcs"] = Int64[a.birth_mcs for a in ar]
    f["archive/end_mcs"] = Int64[a.end_mcs for a in ar]
    f["archive/fate"] = String[string(a.fate) for a in ar]
    f["archive/lifetime_dose_Gy"] = Float64[a.lifetime_dose_Gy for a in ar]
end

"""
    export_transport_snapshot(SR, sim, path; config_toml_path = nothing)

Portable snapshot for the Python/OpenMC side. Carries labels and reference
fields only; asserts NO material composition — that mapping is config-defined
and applied Python-side (fails loudly there when the config is absent).
"""
function export_transport_snapshot(SR, sim, path; config_toml_path = nothing)
    state = sim.state
    species, lineage, generation = _label_arrays(SR, state)
    cfg = config_toml_path === nothing ? "" : read(config_toml_path, String)
    h5open(path, "w") do f
        _common_attrs!(f, state)
        a = attributes(f)
        a["label_state_hash"] =
            label_state_hash(state.lattice, species, lineage, generation)
        a["material_class_source"] = config_toml_path === nothing ?
            "absent" : "config"
        f["config_toml"] = cfg
        f["lattice/cell_id"] = Array{Int32}(state.lattice)
        f["lattice/species_id"] = species
        f["lattice/lineage_id"] = lineage
        f["lattice/generation"] = generation
        f["lattice/interior_mask"] = Array{UInt8}(state.interior)
        f["fields/radiation_cpm"] = state.radiation
        f["fields/melanin"] = state.melanin
        f["dose/accumulated_Gy"] = state.accumulated_dose_Gy
        _write_cells!(f, state)
        _write_events!(f, state)
        f["orientation_probes"] = _orientation_probes(state.lattice)
    end
    return path
end

# Every CPMParams/RadiolysisParams field written individually: these types
# live in the ephemeral load_serial() module, so binary serialization has no
# stable type path across processes. The MersenneTwister (Random stdlib) is
# binary-serialized — stable — and pinned to the writer's Julia version.
function _write_params!(f, group, p)
    for name in fieldnames(typeof(p))
        f["$group/$name"] = getfield(p, name)
    end
end

_read_params(SR, f, group, T) = getfield(SR, T)(;
    (Symbol(k) => read(f["$group/$k"]) for k in keys(f[group]))...)

"""
    export_restart_checkpoint(SR, sim, path)

Complete simulation state for bit-exact resume via `restore_restart_checkpoint`.
"""
function export_restart_checkpoint(SR, sim, path)
    state = sim.state
    h5open(path, "w") do f
        _common_attrs!(f, state)
        a = attributes(f)
        a["sim_mcs"] = sim.mcs
        a["sim_physical_time_s"] = sim.physical_time_s
        a["next_id"] = state.next_id
        a["dose_active"] = Int(state.dose_active)
        a["membrane_dose_rate_Gy_s"] = state.membrane_dose_rate_Gy_s
        a["rd_m"] = sim.rd.m
        a["rd_t"] = sim.rd.t
        f["lattice/cell_id"] = Array{Int32}(state.lattice)
        f["lattice/interior_mask"] = Array{UInt8}(state.interior)
        f["fields/radiation"] = state.radiation
        f["fields/melanin"] = state.melanin
        f["fields/melanin_drive"] = state.melanin_drive
        f["fields/nutrient"] = state.nutrient
        f["dose/rate_mean_Gy_s"] = state.dose_rate_mean_Gy_s
        f["dose/rate_sd_Gy_s"] = state.dose_rate_sd_Gy_s
        f["dose/increment_Gy"] = state.dose_increment_Gy
        f["dose/accumulated_Gy"] = state.accumulated_dose_Gy
        f["contaminant"] = sim.contaminant
        _write_cells!(f, state)
        _write_events!(f, state)
        _write_params!(f, "cpm_params", state.params)
        _write_params!(f, "rd_params", sim.rd.params)
        f["rd/r_grid"] = sim.rd.r_grid
        f["rd/c"] = sim.rd.c
        f["rd/s"] = sim.rd.s
        # RULE 4: the producer declares the semantics.  c and s here are
        # computed on a biomass basis that RADIODIALYSIS: BLOCKED gates
        # whenever X_total != 1.0 -- mean(compute_radial_biomass(...)) is one
        # species' occupied sites over all interior sites, neither a biomass
        # fraction nor a reducer fraction.  A consumer cannot tell that by
        # looking at the arrays, so the file says it.  Asserted by
        # tests/checkpoint_io_tests.jl.
        f["rd/basis_gate_blocked"] = sim.rd.params.X_total != 1.0
        f["rd/basis_gate_note"] = sim.rd.params.X_total != 1.0 ?
            "RADIODIALYSIS: BLOCKED -- c and s were computed on a gated " *
            "biomass basis (X_total = $(sim.rd.params.X_total), a site-occupancy " *
            "mean). No magnitude read from them is a claim about a biofilm." :
            "X_total == 1.0, the standalone default; the basis gate does not apply."
        buf = IOBuffer()
        serialize(buf, sim.rng)
        f["rng/serialized"] = take!(buf)
    end
    return path
end

"""
    restore_restart_checkpoint(SR, path; allow_version_mismatch = false)

Rebuild a `CoupledSimulation` whose continuation is bit-identical to the
unbroken run. Refuses a Julia-version mismatch unless overridden (the RNG
byte stream is version-pinned).
"""
function restore_restart_checkpoint(SR, path; allow_version_mismatch::Bool = false)
    h5open(path, "r") do f
        a = attributes(f)
        ver = read(a["julia_version"])
        if ver != string(VERSION) && !allow_version_mismatch
            error("restart checkpoint written under Julia $ver, running $(VERSION); " *
                  "RNG state is version-pinned (allow_version_mismatch=true to override)")
        end

        params = _read_params(SR, f, "cpm_params", :CPMParams)
        rp = _read_params(SR, f, "rd_params", :RadiolysisParams)

        ids = read(f["cells/id"])
        sp = read(f["cells/species"]); vol = read(f["cells/volume"])
        lin = read(f["cells/lineage_id"]); par = read(f["cells/parent_id"])
        gen = read(f["cells/generation"]); bm = read(f["cells/birth_mcs"])
        ld = read(f["cells/lifetime_dose_Gy"]); com = read(f["cells/com"])
        cells = Dict{Int, Any}()
        for (j, id) in enumerate(ids)
            cells[Int(id)] = SR.CellInfo(Int(sp[j]), Int(vol[j]), com[:, j],
                Int(lin[j]), Int(par[j]), Int(gen[j]), Int(bm[j]), ld[j])
        end
        cells = convert(Dict{Int, SR.CellInfo}, cells)

        dl = [SR.DivisionEvent(
                  Int(read(f["events/division/mcs"])[i]),
                  read(f["events/division/physical_time_s"])[i],
                  Int(read(f["events/division/parent_id"])[i]),
                  Int(read(f["events/division/daughter_a_id"])[i]),
                  Int(read(f["events/division/daughter_b_id"])[i]),
                  Int(read(f["events/division/lineage_id"])[i]),
                  Int(read(f["events/division/parent_generation"])[i]),
                  Int(read(f["events/division/daughter_generation"])[i]),
                  Int(read(f["events/division/parent_volume"])[i]),
                  Int(read(f["events/division/daughter_a_volume"])[i]),
                  Int(read(f["events/division/daughter_b_volume"])[i]))
              for i in eachindex(read(f["events/division/mcs"]))]
        de = [SR.DeathEvent(
                  Int(read(f["events/death/mcs"])[i]),
                  read(f["events/death/physical_time_s"])[i],
                  Int(read(f["events/death/cell_id"])[i]),
                  Int(read(f["events/death/species"])[i]),
                  Int(read(f["events/death/lineage_id"])[i]),
                  Int(read(f["events/death/generation"])[i]),
                  read(f["events/death/lifetime_dose_Gy"])[i])
              for i in eachindex(read(f["events/death/mcs"]))]
        ar = [SR.ArchivedCell(
                  Int(read(f["archive/id"])[i]),
                  Int(read(f["archive/species"])[i]),
                  Int(read(f["archive/lineage_id"])[i]),
                  Int(read(f["archive/parent_id"])[i]),
                  Int(read(f["archive/generation"])[i]),
                  Int(read(f["archive/birth_mcs"])[i]),
                  Int(read(f["archive/end_mcs"])[i]),
                  Symbol(read(f["archive/fate"])[i]),
                  read(f["archive/lifetime_dose_Gy"])[i])
              for i in eachindex(read(f["archive/id"]))]

        state = SR.CPMState(
            Array{Int32}(read(f["lattice/cell_id"])),
            cells,
            read(f["fields/radiation"]),
            read(f["fields/melanin"]),
            read(f["fields/nutrient"]),
            BitArray(read(f["lattice/interior_mask"]) .== 1),
            Int(read(a["next_id"])),
            params,
            SR.build_J_matrix(),
            Int(read(a["mcs"])),
            read(a["physical_time_s"]),
            convert(Vector{SR.DivisionEvent}, dl),
            convert(Vector{SR.DeathEvent}, de),
            convert(Vector{SR.ArchivedCell}, ar),
            read(a["dose_active"]) == 1,
            read(f["dose/rate_mean_Gy_s"]),
            read(f["dose/rate_sd_Gy_s"]),
            read(f["dose/increment_Gy"]),
            read(f["dose/accumulated_Gy"]),
            read(a["membrane_dose_rate_Gy_s"]),
            read(f["fields/melanin_drive"]))

        rd = SR.RadiolysisState(read(f["rd/r_grid"]), read(f["rd/c"]),
                                read(f["rd/s"]), read(a["rd_m"]),
                                read(a["rd_t"]), rp)
        rng = deserialize(IOBuffer(read(f["rng/serialized"])))
        return SR.CoupledSimulation(state, rd, read(f["contaminant"]), rng,
                                    Int(read(a["sim_mcs"])),
                                    read(a["sim_physical_time_s"]))
    end
end

# ---- CLI (demo export from a small fresh run) ----------------

function _cli(args)
    length(args) >= 2 || begin
        println("usage: export_checkpoint.jl <transport|restart> <out.h5> [--config cfg.toml] [--seed N] [--mcs N]")
        exit(1)
    end
    mode, out = args[1], args[2]
    getopt(flag, default) = begin
        i = findfirst(==(flag), args)
        i === nothing ? default : args[i + 1]
    end
    seed = parse(Int, getopt("--seed", "42"))
    n_mcs = parse(Int, getopt("--mcs", "20"))
    cfg = getopt("--config", nothing)

    SR = load_serial()
    Base.invokelatest(_cli_export, SR, mode, out, cfg, seed, n_mcs)
end

function _cli_export(SR, mode, out, cfg, seed, n_mcs)
    sim = SR.init_coupled_simulation(
        SR.CPMParams(N = 20, n_cells_per_species = 2, snapshot_interval = 100),
        # basis_gate_ack: this CLI exports an interchange artifact rather than
        # asserting anything, and every file it writes carries
        # rd/basis_gate_blocked so a consumer cannot read c or s without seeing
        # that the basis was gated. Enumerated in the ack census in
        # tests/radiodialysis_basis_gate.jl.
        SR.RadiolysisParams(Nr = 20, Ddot_R = 1.0, c_ext = 1.0,
                            basis_gate_ack = true); seed)
    SR.advance_window!(sim, n_mcs)
    if mode == "transport"
        export_transport_snapshot(SR, sim, out; config_toml_path = cfg)
    elseif mode == "restart"
        export_restart_checkpoint(SR, sim, out)
    else
        error("unknown mode $mode")
    end
    @printf("wrote %s (%s, mcs=%d, seed=%d)\n", out, mode, n_mcs, seed)
end

if abspath(PROGRAM_FILE) == @__FILE__
    _cli(ARGS)
end
