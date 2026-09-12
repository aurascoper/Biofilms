# The renderer-free half of the viewer: reading a transport snapshot into the UInt8 grid
# GLMakie's voxels recipe draws, the palette and labels shared with the serial script's
# figures, and the CLI option parser. Needs HDF5 only, so tests/vti_export_tests.jl includes
# it from the main environment and exercises it in CI; viewer/visualize_lattice.jl includes it
# and adds the OpenGL half.
using HDF5

const COLORS = ["#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6"]
const LABELS = ["C. neoformans", "D. radiodurans", "C. sphaerospermum",
                "B. subtilis", "A. niger", "S. oneidensis", "O. intermedium"]

"""
    species_grid(snapshot) -> (UInt8 array, mcs)

0x00 where the site is medium or wall (cell_id equal to the file's own
`cell_id_background` or `cell_id_wall`), the species index 1..7 elsewhere. A file whose
`logical_axis_order` is not "xyz" is refused: this grid maps Julia axis 1 to x.
"""
function species_grid(snapshot::AbstractString)
    h5open(snapshot, "r") do f
        a = HDF5.attributes(f)
        order = read(a["logical_axis_order"])
        order == "xyz" || throw(ArgumentError("logical_axis_order=\"$order\"; this viewer maps axis 1 to x and knows only \"xyz\""))
        background = read(a["cell_id_background"]); wall = read(a["cell_id_wall"])
        cell_id = read(f["lattice/cell_id"]); species_id = read(f["lattice/species_id"])
        grid = map((c, s) -> (c == background || c == wall) ? 0x00 : UInt8(s), cell_id, species_id)
        return grid, Int(read(a["mcs"]))
    end
end

"""
    viewer_options(args) -> (snapshot, still, record_to, frames)

The CLI contract: the first argument is the snapshot; `--still out.png` or `--record out.mp4`
choose an output, `--frames N` the orbit length. A flag without a value, an unknown flag, or
both outputs at once is refused.
"""
function viewer_options(args::AbstractVector{<:AbstractString})
    isempty(args) && throw(ArgumentError("usage: visualize_lattice.jl <transport_snapshot.h5> [--still out.png | --record out.mp4] [--frames N]"))
    snapshot = args[1]
    opts = Dict{String, String}()
    i = 2
    while i <= length(args)
        flag = args[i]
        flag in ("--still", "--record", "--frames") || throw(ArgumentError("unknown option $flag"))
        i + 1 <= length(args) || throw(ArgumentError("$flag needs a value"))
        opts[flag] = args[i + 1]
        i += 2
    end
    haskey(opts, "--still") && haskey(opts, "--record") && throw(ArgumentError("--still and --record are exclusive"))
    frames = parse(Int, get(opts, "--frames", "120"))
    frames >= 1 || throw(ArgumentError("--frames must be >= 1"))
    return (snapshot = snapshot, still = get(opts, "--still", nothing),
            record_to = get(opts, "--record", nothing), frames = frames)
end
