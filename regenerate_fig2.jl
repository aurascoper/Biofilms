#!/usr/bin/env julia
# Regenerate fig2_melanin_accumulation ONLY.
#
#   julia --project=. regenerate_fig2.jl
#   julia --project=. regenerate_fig2.jl --self-test    # the refusal must bite
#
# WHY THIS EXISTS. The committed figure printed an ordering attribution inside the image,
# which SEED42-SCOPE-01 requalified on 2026-09-17: the ordering fails in 19 of 256 seeds,
# so the input does not fix it. The row holds the wording. A caption inside an image is the
# sink a prose audit cannot reach (FIG-01 to FIG-07 record that failure), so the
# generator's text changed and the artifact has to be rebuilt to match it.
#
# WHY NOT `julia --project=. biofilms_potts.jl`. Its entry point is main_coupled(), which
# cannot run while RADIODIALYSIS: BLOCKED stands. export_figures writes figs 1 and 2
# together, so this renders into a temporary directory and copies exactly one figure
# across, as regenerate_fig3.jl does.
#
# WHAT THE FIGURE PLOTS IS PINNED, EVERY POINT OF IT. The committed reference is
# diagnostics/melanin_ensemble/sweep_n40_p6_seed42_int20.csv, written by this repository's
# own sweep.jl at the figure's configuration and snapshot interval. Every plotted snapshot
# of every species must match it, and the species set must be complete: a truncated
# reference that happens to hold one row would otherwise validate one number and pass the
# rest. The 16-seed sweep is read as well, so the two committed artifacts must agree at
# MCS 100. The script REFUSES to copy on any mismatch, and --self-test proves the refusal
# bites by mutating a value taken from the reference itself.
using Printf

const HERE   = @__DIR__
const OUT    = joinpath(HERE, "preprint", "figures")
const NAME   = "fig2_melanin_accumulation"
const ENS    = joinpath(HERE, "diagnostics", "melanin_ensemble")
const REF    = joinpath(ENS, "sweep_n40_p6_seed42_int20.csv")
const SWEEP  = joinpath(ENS, "sweep_n40_p6_seeds42-57.csv")
const MCS    = 100
const STEP   = 20
const SEED   = 42

M = Module(:Fig2Regen)
Base.eval(M, :(using LinearAlgebra, Statistics, Random, Printf))
Base.include(M, joinpath(HERE, "biofilms_potts.jl"))

const CPMParams      = Base.eval(M, :CPMParams)
const SPECIES_NAMES  = Base.eval(M, :SPECIES_NAMES)
const run_simulation = Base.eval(M, :run_simulation)
const export_figures = Base.eval(M, :export_figures)

"""(mcs, species name) => mean melanin, for seed SEED, from a committed sweep CSV."""
function reference(path::String)
    out = Dict{Tuple{Int,String},Float64}()
    for line in eachline(path)
        (startswith(line, "#") || startswith(line, "seed,")) && continue
        f = split(line, ",")
        parse(Int, f[1]) == SEED || continue
        out[(parse(Int, f[2]), String(f[4]))] = parse(Float64, f[8])
    end
    return out
end

"""Every way the reference could be incomplete, named rather than assumed."""
function reference_gaps(ref::Dict{Tuple{Int,String},Float64})
    gaps = String[]
    for mcs in 0:STEP:MCS, name in SPECIES_NAMES
        haskey(ref, (mcs, name)) || push!(gaps, "MCS $mcs, $name")
    end
    extra = setdiff(Set(k[2] for k in keys(ref)), Set(SPECIES_NAMES))
    isempty(extra) || push!(gaps, "names the model does not have: $(sort(collect(extra)))")
    return gaps
end

"""Every plotted point that differs from the reference."""
function mismatches(traj, ref::Dict{Tuple{Int,String},Float64})
    bad = String[]
    plotted = [s for s in traj if s.mcs % STEP == 0 && s.mcs <= MCS]
    [s.mcs for s in plotted] == collect(0:STEP:MCS) ||
        push!(bad, "the plotted snapshots are $([s.mcs for s in plotted]), not $(collect(0:STEP:MCS))")
    for s in plotted, sd in s.species_data
        want = ref[(s.mcs, sd.name)]
        abs(sd.mean_melanin - want) <= 5e-7 ||
            push!(bad, @sprintf("MCS %d, %s: rendered %.6f, committed %.6f", s.mcs, sd.name, sd.mean_melanin, want))
    end
    return bad
end

function render()
    params = CPMParams(N = 40, n_cells_per_species = 6, snapshot_interval = STEP)
    _, traj = run_simulation(params, MCS; seed = SEED)
    return params, traj
end

ref = reference(REF)
gaps = reference_gaps(ref)
isempty(gaps) || error("REFUSING: the reference in $REF does not cover $(length(gaps)) of the " *
                       "plotted points: $(first(gaps, 5))")

if "--self-test" in ARGS
    # The known-bad is taken from the artifact, never typed: one value from the reference,
    # moved by 1e-3. A mutation that lands nowhere would pass silently, so the mutated
    # reference is asserted to differ first.
    params, traj = render()
    isempty(mismatches(traj, ref)) || error("SELF-TEST: the unmutated reference already mismatches")
    victim = (STEP, SPECIES_NAMES[3])
    mutated = copy(ref)
    mutated[victim] = ref[victim] + 1e-3
    mutated[victim] == ref[victim] && error("SELF-TEST: the mutation did not apply at $victim")
    bad = mismatches(traj, mutated)
    length(bad) == 1 && occursin(SPECIES_NAMES[3], bad[1]) ||
        error("SELF-TEST: a mutated reference produced $(length(bad)) mismatches: $bad")
    println("  self-test: a mutated reference value is refused, and names $(SPECIES_NAMES[3]) at MCS $STEP")
    println("  self-test: the unmutated reference passes, so the check is not refusing everything")
    exit(0)
end

# The two committed sweeps must agree where they overlap, at MCS 100.
overlap = reference(SWEEP)
disagree = [n for n in SPECIES_NAMES if abs(get(overlap, (MCS, n), NaN) - ref[(MCS, n)]) > 5e-7]
isempty(disagree) || error("REFUSING: the two committed sweeps disagree at MCS $MCS for $disagree")

params, traj = render()
bad = mismatches(traj, ref)
for mcs in 0:STEP:MCS, n in SPECIES_NAMES[[1, 3, 5]]
    @printf("  MCS %3d  %-20s committed %.6f\n", mcs, n, ref[(mcs, n)])
end
isempty(bad) || error("REFUSING TO COPY: the render does not match the committed reference:\n  " *
                      join(bad, "\n  ") * "\nThe figure would plot numbers no committed run supports.")

tmp = mktempdir(; cleanup = false)
export_figures(traj, nothing, params; outdir = tmp)
for ext in ("pdf", "png")
    cp(joinpath(tmp, "$NAME.$ext"), joinpath(OUT, "$NAME.$ext"); force = true)
end
rm(tmp; recursive = true)

@printf("\n  every plotted point matches the committed reference; copied %s.{pdf,png} into %s\n", NAME, OUT)
println("\n  Now refresh the sidecars:")
println("    cd preprint/figures")
println("    pdftotext -layout $NAME.pdf $NAME.txt")
println("    sha256sum $NAME.pdf | cut -d' ' -f1 > $NAME.sha256")
