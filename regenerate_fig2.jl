#!/usr/bin/env julia
# Regenerate fig2_melanin_accumulation ONLY.
#
#   julia --project=. regenerate_fig2.jl
#
# WHY THIS EXISTS. The committed figure printed "Ordering follows the
# hand-specified alpha_M_species scales" inside the image. SEED42-SCOPE-01
# requalified that on 2026-09-17: the ordering fails in 19 of 256 seeds, so the
# input does not fix it. A caption inside an image is the sink a prose audit
# cannot reach (FIG-01 to FIG-07 record the same failure), so the generator's
# text changed and the artifact has to be rebuilt to match it.
#
# WHY NOT `julia --project=. biofilms_potts.jl`. Its entry point is
# main_coupled(), which cannot run while RADIODIALYSIS: BLOCKED stands.
# export_figures writes figs 1 and 2 together, so this renders into a temporary
# directory and copies exactly one figure across, as regenerate_fig3.jl does.
#
# THE CONFIGURATION IS NOT A GUESS. The committed figure ends at 1.44 for
# C. sphaerospermum, and diagnostics/melanin_ensemble/sweep_n40_p6_seeds42-57.csv
# gives 1.437157 for seed 42 at MCS 100 under N=40, 6 parcels per species. This
# script runs that configuration and REFUSES to copy unless all three producers
# match that committed ensemble at MCS 100. A re-render that quietly plots other
# numbers is the failure this check exists to stop.

using Printf

const HERE   = @__DIR__
const OUT    = joinpath(HERE, "preprint", "figures")
const NAME   = "fig2_melanin_accumulation"
const SWEEP  = joinpath(HERE, "diagnostics", "melanin_ensemble", "sweep_n40_p6_seeds42-57.csv")
const MCS    = 100
const SEED   = 42

M = Module(:Fig2Regen)
Base.eval(M, :(using LinearAlgebra, Statistics, Random, Printf))
Base.include(M, joinpath(HERE, "biofilms_potts.jl"))

const CPMParams      = Base.eval(M, :CPMParams)
const run_simulation = Base.eval(M, :run_simulation)
const export_figures = Base.eval(M, :export_figures)

"""Committed mean melanin at MCS 100, per species name, for seed SEED."""
function committed_at_mcs()
    out = Dict{String,Float64}()
    for line in eachline(SWEEP)
        startswith(line, "#") && continue
        f = split(line, ",")
        f[1] == "seed" && continue
        parse(Int, f[1]) == SEED && parse(Int, f[2]) == MCS || continue
        out[f[4]] = parse(Float64, f[8])
    end
    isempty(out) && error("no seed $SEED rows at MCS $MCS in $SWEEP")   # plumbing, not the check
    return out
end

params = CPMParams(N = 40, n_cells_per_species = 6, snapshot_interval = 20)
_, traj = run_simulation(params, MCS; seed = SEED)

final = last(traj)
final.mcs == MCS || error("the last snapshot is MCS $(final.mcs), not $MCS")
rendered = Dict(sd.name => sd.mean_melanin for sd in final.species_data)
committed = committed_at_mcs()

bad = String[]
for (name, value) in sort(collect(committed))
    got = get(rendered, name, NaN)
    @printf("  %-24s rendered %.6f   committed %.6f\n", name, got, value)
    abs(got - value) <= 5e-7 || push!(bad, name)
end
isempty(bad) || error("REFUSING TO COPY: these species do not match the committed ensemble at " *
                      "MCS $MCS: $bad. The figure would plot numbers no committed run supports.")

tmp = mktempdir(; cleanup = false)
export_figures(traj, nothing, params; outdir = tmp)
for ext in ("pdf", "png")
    cp(joinpath(tmp, "$NAME.$ext"), joinpath(OUT, "$NAME.$ext"); force = true)
end
rm(tmp; recursive = true)

@printf("\n  copied %s.{pdf,png} into %s\n", NAME, OUT)
println("\n  Now refresh the sidecars:")
println("    cd preprint/figures")
println("    pdftotext -layout $NAME.pdf $NAME.txt")
println("    sha256sum $NAME.pdf | cut -d' ' -f1 > $NAME.sha256")
