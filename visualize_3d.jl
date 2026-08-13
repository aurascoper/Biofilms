#!/usr/bin/env julia
# 3D voxel animation of the serial Potts biofilm (biofilms_potts.jl).
# Cubes = occupied lattice sites, colored by species; camera orbits as
# the biofilm evolves. Output: biofilm_3d.gif in the repo root.
#
#   julia --project=. visualize_3d.jl [n_mcs] [interval] [seed]
#   defaults: 100 5 42

using Random, Printf

n_mcs    = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 100
interval = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 5
seed     = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 42

# Same loading trick as validate_serial.jl: serial script minus its figure section.
function load_serial()
    src = read(joinpath(@__DIR__, "biofilms_potts.jl"), String)
    src = split(src, "#  13. Figure export")[1]
    SerialRef = Module(:SerialRef)
    Base.eval(SerialRef, :(using LinearAlgebra, Statistics, Random, Printf))
    Base.include_string(SerialRef, src, "biofilms_potts.jl")
    return SerialRef
end

# Run the sim (same per-MCS calls as run_simulation), capturing a
# species-index grid (0 = medium/wall) every `interval` MCS.
function simulate(SR::Module, n_mcs::Int, interval::Int, seed::Int)
    params = SR.CPMParams(N = 40, n_cells_per_species = 6)
    rng = MersenneTwister(seed)
    state = SR.init_state(params; seed)
    N = params.N

    spgrid() = UInt8[SR.species_of(state.lattice[x, y, z], state.cells)
                     for x in 1:N, y in 1:N, z in 1:N]

    frames = [(0, spgrid())]
    for mcs in 1:n_mcs
        SR.mcs_step!(state, rng)
        SR.update_melanin!(state)
        SR.update_nutrient!(state)
        if mcs % interval == 0 || mcs == n_mcs
            push!(frames, (mcs, spgrid()))
            @printf("# mcs %d captured\n", mcs)
        end
    end
    return frames
end

frames = Base.invokelatest(simulate, load_serial(), n_mcs, interval, seed)

using CairoMakie

# Palette + labels: same as serial section 13 (FIG_COLORS / FIG_LABELS)
const COLORS = ["#e6194b", "#3cb44b", "#4363d8", "#f58231",
                "#911eb4", "#42d4f4", "#f032e6"]
const LABELS = ["C. neoformans", "D. radiodurans", "C. sphaerospermum",
                "B. subtilis", "A. niger", "S. oneidensis", "O. intermedium"]

N = size(frames[1][2], 1)
fig = Figure(size = (900, 800))
ax = Axis3(fig[1, 1], aspect = :data, azimuth = 0.4, elevation = 0.5,
           limits = (1, N, 1, N, 1, N))
Legend(fig[1, 2],
       [MarkerElement(marker = :rect, markersize = 14, color = c) for c in COLORS],
       LABELS, framevisible = false)
cube = Rect3f(Vec3f(-0.5), Vec3f(1))

record(fig, joinpath(@__DIR__, "biofilm_3d.gif"), enumerate(frames);
       framerate = 6) do (i, (mcs, grid))
    empty!(ax)
    for s in 1:7
        pts = Point3f.(Tuple.(findall(==(UInt8(s)), grid)))
        isempty(pts) && continue
        meshscatter!(ax, pts; marker = cube, markersize = 1,
                     color = COLORS[s], shading = NoShading)
    end
    ax.title = @sprintf("Radiotrophic biofilm — MCS %d", mcs)
    ax.azimuth = 0.4 + 2π * (i - 1) / length(frames)
end
println("wrote biofilm_3d.gif")
