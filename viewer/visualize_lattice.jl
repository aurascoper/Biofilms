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
# docs/visualization/lattice_viewer.md rather than hidden behind a skip. What can be
# checked without a context (the grid, the palette, the CLI contract, and that every
# name this file uses resolves) lives in viewer/lattice_grid.jl and is tested in CI.
#
# Palette and labels are the serial script's FIG_COLORS / FIG_LABELS
# (biofilms_potts.jl, section 13), so species s has the same colour here as in
# every figure the repository ships.
using GLMakie
include(joinpath(@__DIR__, "lattice_grid.jl"))

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
    o = viewer_options(ARGS)
    out = show_lattice(o.snapshot; record_to = o.record_to, still = o.still, frames = o.frames)
    out isa String ? println("wrote $out") : (println("close the window to exit"); wait(GLMakie.Screen()))
end
