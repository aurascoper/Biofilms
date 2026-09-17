#!/usr/bin/env julia
# Regenerate fig3_membrane_transport ONLY.
#
#   julia --project=. regenerate_fig3.jl
#
# WHY THIS EXISTS RATHER THAN `julia --project=. biofilms_potts.jl`. The script's
# default entry point is main_coupled(), and it cannot run: the coupled loop
# reconstructs X_total from occupancy and RADIODIALYSIS: BLOCKED refuses to
# integrate. So figs 3-4 have not been regenerable since the gate landed, and
# 1b51126 -- which regenerated figs 1-2 through main()'s --no-radiolysis path --
# could only add .sha256/.txt sidecars for figs 3-4 over pdfs it could not
# rebuild. A generator nobody can run is the same liability that commit was
# written about, so this is committed rather than kept in a scratch directory.
#
# WHY FIG3 AND NOT FIG4. export_figures writes both, so this renders into a
# temporary directory and copies exactly one file across.
#
#   fig3 plots m(t) and P_eff/P0.  dm/dt = -k_dam*Ddot_R*m and
#   P_eff = P0*exp(alpha_P*Ddot_R*t).  Neither reads X_total or X_red.
#   fig4 plots c_wall, c_mean, s_mean.  The gated basis enters through
#   uptake = k_ads*X_total + k_red*X_red, which drives c and s and nothing else.
#
# basis_gate_ack: the exemption is that claim, and NOTHING HERE TAKES IT ON
# TRUST. The run is repeated at 10x the uptake constants and fig3's two series
# must come back identical while fig4's must NOT -- a control that fails in both
# directions, so it cannot pass by measuring nothing. The copy is refused if
# either half is violated. Mirrors the "exemption's claim is true" testset in
# tests/radiodialysis_basis_gate.jl.

using Printf

const HERE = @__DIR__
const OUT  = joinpath(HERE, "preprint", "figures")
const NAME = "fig3_membrane_transport"

M = Module(:Fig3Regen)
Base.eval(M, :(using LinearAlgebra, Statistics, Random, Printf))
Base.include(M, joinpath(HERE, "biofilms_potts.jl"))

const CPMParams              = Base.eval(M, :CPMParams)
const RadiolysisParams       = Base.eval(M, :RadiolysisParams)
const run_simulation_coupled = Base.eval(M, :run_simulation_coupled)
const export_figures         = Base.eval(M, :export_figures)

function render(outdir; k_ads = 0.05, k_red = 0.02)
    params = CPMParams(N = 40, n_cells_per_species = 6, snapshot_interval = 20)
    rp = RadiolysisParams(Nr = 40, Ddot_R = 1.0, c_ext = 1.0,
                          k_ads = k_ads, k_red = k_red, basis_gate_ack = true)
    _, _, _, ctraj = run_simulation_coupled(params, rp, 100; seed = 42)
    export_figures(ctraj, ctraj, params; outdir)
    return ctraj
end

base_dir = mktempdir(; cleanup = false)
pert_dir = mktempdir(; cleanup = false)
base = render(base_dir)
pert = render(pert_dir; k_ads = 0.5, k_red = 0.2)

series(ct, f) = [f(cs) for cs in ct]
fig3_invariant = series(base, cs -> cs.m)      == series(pert, cs -> cs.m) &&
                 series(base, cs -> cs.P_eff)  == series(pert, cs -> cs.P_eff)
fig4_moves     = series(base, cs -> cs.c_wall) != series(pert, cs -> cs.c_wall)

@printf("\n  fig3 series invariant to 10x uptake : %s  (required true)\n", fig3_invariant)
@printf("  fig4 series moves under 10x uptake  : %s  (required true -- the control must bite)\n",
        fig4_moves)

fig3_invariant || error("REFUSING TO COPY: m or P_eff moved with the gated basis, so fig3 " *
                        "does record a basis-dependent quantity and the exemption is void.")
fig4_moves     || error("REFUSING TO COPY: the control did not bite -- c_wall was unchanged " *
                        "by a 10x uptake perturbation, so it proves nothing about fig3.")

for ext in ("pdf", "png")
    cp(joinpath(base_dir, "$NAME.$ext"), joinpath(OUT, "$NAME.$ext"); force = true)
end
rm(base_dir; recursive = true); rm(pert_dir; recursive = true)

@printf("\n  copied %s.{pdf,png} into %s\n", NAME, OUT)
@printf("  m(t)     = %s\n", round.(series(base, cs -> cs.m), digits = 4))
@printf("  P_eff/P0 = %s\n",
        round.(series(base, cs -> cs.P_eff) ./ series(base, cs -> cs.P_eff)[1], digits = 4))
println("\n  Now refresh the sidecars:")
println("    cd preprint/figures")
println("    pdftotext -layout $NAME.pdf $NAME.txt")
println("    sha256sum $NAME.pdf | cut -d' ' -f1 > $NAME.sha256")
