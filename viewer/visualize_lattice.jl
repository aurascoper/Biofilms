#!/usr/bin/env julia
# Interactive 3-D view of a transport snapshot: one cube per occupied lattice
# site, coloured by species, on GLMakie's `voxels` recipe. Reads the snapshot;
# never re-runs the simulation. Lattice units throughout: the axes are site
# indices, the title carries the MCS and nothing else.
#
#   julia --project=viewer viewer/visualize_lattice.jl <transport_snapshot.h5> [--still out.png | --record out.mp4] [--frames N]
#
# This project pins GLMakie separately (viewer/Project.toml, viewer/Manifest.toml)
# so the main environment pulls no OpenGL dependency. The viewer is not run in CI:
# GLMakie needs an OpenGL 3.3 context, and that is stated as uncovered surface in
# docs/visualization/lattice_viewer.md rather than hidden behind a skip.
#
# Palette and labels are the serial script's FIG_COLORS / FIG_LABELS
# (biofilms_potts.jl, section 13), so species s has the same colour here as in
# every figure the repository ships.
using HDF5, GLMakie

const COLORS = ["#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6"]
const LABELS = ["C. neoformans", "D. radiodurans", "C. sphaerospermum",
                "B. subtilis", "A. niger", "S. oneidensis", "O. intermedium"]

"""
    species_grid(snapshot) -> (UInt8 array, mcs)

0x00 where the site is medium or wall (cell_id equal to the file's own
`cell_id_background` or `cell_id_wall`), the species index 1..7 elsewhere.
"""
function species_grid(snapshot::AbstractString)
    h5open(snapshot, "r") do f
        a = HDF5.attributes(f)
        read(a["logical_axis_order"]) == "xyz" || error("logical_axis_order is not xyz; this viewer maps axis 1 to x")
        background = read(a["cell_id_background"]); wall = read(a["cell_id_wall"])
        cell_id = read(f["lattice/cell_id"]); species_id = read(f["lattice/species_id"])
        grid = map((c, s) -> (c == background || c == wall) ? 0x00 : UInt8(s), cell_id, species_id)
        return grid, Int(read(a["mcs"]))
    end
end

function show_lattice(snapshot::AbstractString; record_to = nothing, still = nothing, frames::Int = 120)
    grid, mcs = species_grid(snapshot)
    N = size(grid)
    present = sort(unique(filter(!=(0x00), grid)))
    fig = Figure(size = (1000, 800))
    ax = Axis3(fig[1, 1]; aspect = :data, title = "MCS $mcs",
               xlabel = "x (sites)", ylabel = "y (sites)", zlabel = "z (sites)",
               limits = (0, N[1], 0, N[2], 0, N[3]))
    # voxels: value 0x00 is air and is never drawn; ids 1..7 index `color`.
    voxels!(ax, 0 .. N[1], 0 .. N[2], 0 .. N[3], grid;
            color = parse.(Makie.Colorant, COLORS), is_air = ==(0x00))
    Legend(fig[1, 2],
           [PolyElement(color = COLORS[s]) for s in present],
           [LABELS[s] for s in present]; framevisible = false)
    if still !== nothing
        save(still, fig)
        return still
    end
    if record_to === nothing
        display(fig)
        return fig
    end
    record(fig, record_to, 1:frames; framerate = 24) do i
        ax.azimuth[] = 2π * (i - 1) / frames
    end
    return record_to
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && (println("usage: visualize_lattice.jl <transport_snapshot.h5> [--still out.png | --record out.mp4] [--frames N]"); exit(1))
    opt(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i + 1])
    out = show_lattice(ARGS[1]; record_to = opt("--record", nothing), still = opt("--still", nothing),
                       frames = parse(Int, opt("--frames", "120")))
    out isa String ? println("wrote $out") : (println("close the window to exit"); wait(GLMakie.Screen()))
end
