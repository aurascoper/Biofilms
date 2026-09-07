#!/usr/bin/env julia
# transport_snapshot.h5 -> .vti (VTK ImageData) for ParaView, and a directory
# of snapshots -> .pvd keyed by each snapshot's `mcs` attribute.
#
# Lattice units only. The spacing written is 1.0 per site unless a declared
# pitch is passed explicitly, because the physical pitch is D-PITCH and the
# clock is D-TIMESERIES, both `awaiting_measurement` in
# data/calibration/reference_d_requirements.csv. The snapshot's own axis and
# sentinel declarations (docs/exchange_schema.md) are read from its attributes
# and carried into the file as field data, never restated here.
#
#   julia --project=. export_vti.jl <snapshot.h5> <out_stem> [--restart r.h5] [--dose d.h5]
#   julia --project=. export_vti.jl <snapshot_dir> <out_stem> [--restart-dir d] [--dose-dir d]
#
# Cell data (one value per lattice site, VTK x = Julia axis 1):
#   species          UInt8   0 = medium or wall, 1..7 = species index
#   cell_id          Int32   as stored: background and wall sentinels kept
#   lineage_id, generation, interior_mask, radiation_cpm, melanin  as stored
#   accumulated_dose_Gy      as stored; physical Gy, zero until a dose was imported
#   nutrient                 only with --restart (the snapshot does not carry it)
#   dose_rate_mean_Gy_s      only with --dose (transport_result_*.h5, mesh/dose_rate_mean_Gy_s), labelled
#                            from that file's own source_rate_photons_per_s and target_calibration attributes
using HDF5, WriteVTK

const CARRIED_NUMERIC = ("schema_version", "coordinate_index_base", "cell_id_background",
                         "cell_id_wall", "mcs", "physical_time_s")
const CARRIED_STRING = ("logical_axis_order", "dataset_axis_order_h5py", "git_sha")
const SITE_ARRAYS = (("lattice/cell_id", "cell_id"), ("lattice/lineage_id", "lineage_id"),
                     ("lattice/generation", "generation"), ("lattice/interior_mask", "interior_mask"),
                     ("fields/radiation_cpm", "radiation_cpm"), ("fields/melanin", "melanin"),
                     ("dose/accumulated_Gy", "accumulated_dose_Gy"))

"""
    export_vti(snapshot, stem; restart=nothing, dose=nothing, spacing=1.0, declared_pitch=nothing)

Write `stem.vti` from a transport snapshot. `spacing != 1.0` is refused unless
`declared_pitch = (value=..., unit=..., source=...)` is passed, in which case the
spacing is the declared value and the file says so in its `units` field.
Returns the WriteVTK dataset (saved unless `save=false`).
"""
function export_vti(snapshot::AbstractString, stem::AbstractString;
                    restart = nothing, dose = nothing, spacing::Real = 1.0,
                    declared_pitch = nothing, save::Bool = true)
    if spacing != 1.0 && declared_pitch === nothing
        throw(ArgumentError("spacing=$spacing is a physical length per site, and the lattice pitch " *
              "is D-PITCH, awaiting_measurement (data/calibration/reference_d_requirements.csv). " *
              "Pass declared_pitch=(value=..., unit=..., source=...) to write a declared spacing."))
    end
    if declared_pitch !== nothing
        spacing = declared_pitch.value
        units = "declared: $(declared_pitch.value) $(declared_pitch.unit) per site; source: $(declared_pitch.source)"
    else
        units = "lattice"
    end
    h5open(snapshot, "r") do f
        a = attributes(f)
        order = read(a["logical_axis_order"])
        order == "xyz" || throw(ArgumentError("logical_axis_order=\"$order\": this exporter maps Julia " *
                                              "axis 1 to VTK x and knows only \"xyz\""))
        cell_id = read(f["lattice/cell_id"])
        N = size(cell_id)
        background = read(a["cell_id_background"]); wall = read(a["cell_id_wall"])
        species_id = read(f["lattice/species_id"])
        species = map((c, s) -> (c == background || c == wall) ? 0x00 : UInt8(s), cell_id, species_id)
        axes = ntuple(i -> range(0.0, step = Float64(spacing), length = N[i] + 1), 3)
        vtk = vtk_grid(stem, axes...; compress = false)
        vtk["species", VTKCellData()] = species
        for (ds, name) in SITE_ARRAYS
            vtk[name, VTKCellData()] = read(f[ds])
        end
        if restart !== nothing
            nut = h5open(g -> read(g["fields/nutrient"]), restart, "r")
            size(nut) == N || throw(ArgumentError("restart nutrient field $(size(nut)) is not the lattice $N"))
            vtk["nutrient", VTKCellData()] = nut
        end
        if dose !== nothing
            d, label = h5open(dose, "r") do g
                da = attributes(g)
                # The file's own qualifiers, or nothing: a Gy/s array without them is the
                # withdrawn-annotation class with a different unit.
                haskey(da, "target_calibration") || throw(ArgumentError(
                    "$dose carries no target_calibration attribute; a Gy/s field cannot be labelled without it"))
                haskey(da, "source_rate_photons_per_s") || throw(ArgumentError(
                    "$dose carries no source_rate_photons_per_s attribute"))
                target = Bool(read(da["target_calibration"]))
                rate = read(da["source_rate_photons_per_s"])
                qualifier = target ? "target_calibration = true" :
                    "target_calibration = false: synthetic source rate, not a physical target"
                read(g["mesh/dose_rate_mean_Gy_s"]),
                "Gy s^-1, schema mesh/dose_rate_mean_Gy_s, at source_rate_photons_per_s = $rate; $qualifier"
            end
            size(d) == N || throw(ArgumentError("dose mesh $(size(d)) is not the lattice $N; this exporter " *
                                                "resamples nothing (viewer_bundle.h5 is where that happens)"))
            vtk["dose_rate_mean_Gy_s", VTKCellData()] = d
            vtk["dose_rate_mean_Gy_s_units", VTKFieldData()] = label
        end
        vtk["units", VTKFieldData()] = units
        vtk["species_zero", VTKFieldData()] = "0 = medium or wall (cell_id $background or $wall); 1..7 = species index"
        vtk["accumulated_dose_Gy_units", VTKFieldData()] =
            "Gy, physical, schema dose/accumulated_Gy: zeros until a dose was imported"
        for k in CARRIED_NUMERIC
            vtk[k, VTKFieldData()] = Float64(read(a[k]))
        end
        for k in CARRIED_STRING
            vtk[k, VTKFieldData()] = String(read(a[k]))
        end
        save && vtk_save(vtk)
        return vtk
    end
end

"""
    export_series(dir, stem; restart_dir=nothing, dose_dir=nothing, kwargs...)

Every `*.h5` in `dir` becomes `stem_mcsNNNNNN.vti`, collected in `stem.pvd` with the
snapshot's `mcs` attribute as the time key (a step count, not seconds). Two snapshots
with one `mcs` are refused. Companion files are matched by basename in the other dirs.
"""
function export_series(dir::AbstractString, stem::AbstractString;
                       restart_dir = nothing, dose_dir = nothing, kwargs...)
    files = sort(filter(p -> endswith(p, ".h5"), readdir(dir; join = true)))
    isempty(files) && throw(ArgumentError("no .h5 snapshots in $dir"))
    seen = Dict{Int, String}()
    written = String[]
    paraview_collection(stem) do pvd
        for p in files
            m = h5open(g -> Int(read(attributes(g)["mcs"])), p, "r")
            haskey(seen, m) && throw(ArgumentError("mcs=$m is carried by both $(seen[m]) and $p"))
            seen[m] = p
            companion(d) = d === nothing ? nothing : joinpath(d, basename(p))
            vtk = export_vti(p, "$(stem)_mcs$(lpad(m, 6, '0'))";
                             restart = companion(restart_dir), dose = companion(dose_dir), kwargs...)
            pvd[Float64(m)] = vtk
            push!(written, "$(stem)_mcs$(lpad(m, 6, '0')).vti")
        end
    end
    return written
end

function _cli(args)
    length(args) >= 2 || begin
        println("usage: export_vti.jl <snapshot.h5|snapshot_dir> <out_stem> [--restart r.h5|--restart-dir d] [--dose d.h5|--dose-dir d]")
        exit(1)
    end
    src, stem = args[1], args[2]
    opt(flag) = (i = findfirst(==(flag), args); i === nothing ? nothing : args[i + 1])
    if isdir(src)
        w = export_series(src, stem; restart_dir = opt("--restart-dir"), dose_dir = opt("--dose-dir"))
        println("wrote $(length(w)) .vti + $stem.pvd (time key = mcs, lattice units)")
    else
        export_vti(src, stem; restart = opt("--restart"), dose = opt("--dose"))
        println("wrote $stem.vti (lattice units; spacing 1.0 per site)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    _cli(ARGS)
end
