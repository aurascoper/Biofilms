#!/usr/bin/env julia
# Fig 5 — acceptance regime of the JACC checkerboard CPM.
#
#   julia --project=. jacc_acceptance_figure.jl
#
# TWO DISTRIBUTIONS, NO MAP. d404438 refused a spatial map for the per-voxel
# decisive-label tally -- at 400 MCS only 6157 of 64000 voxels received any
# accepted move, median 3 -- and the same sparsity argument applies to anything
# drawn per voxel here. What the lattice picture cannot answer is whether the
# Metropolis is in a usable regime at all: acceptance is 1 for ΔH <= 0 and
# exp(-ΔH/T) otherwise, so T_cpm sets everything, and a rate pinned near 0 or
# near 1 renders as a simulation that is plainly running. These two panels are
# that diagnostic.
#
# basis_gate_ack: this script records NO radiodialysis quantity -- it plots
# acceptance rates and a ΔH distribution. `delta_H` reads
# lat/vols/spec/J/beta/melc/rad/mel and never `nut`, so the gated basis cannot
# reach an accepted move; tests/jacc_parity_tests.jl holds that claim to account
# by requiring 10x the uptake constants to leave the tables byte-identical.

using CairoMakie, Statistics, Printf

const HERE = @__DIR__
const OUT  = joinpath(HERE, "preprint", "figures")

Port = Module(:FigPort)
Base.eval(Port, :(using LinearAlgebra, Statistics, Random, Printf))
Base.include(Port, joinpath(HERE, "biofilms_potts_jacc.jl"))

const RUN = Base.eval(Port, :run_coupled)
const RP  = Base.eval(Port, :RadiolysisParams)

# Provenance, printed into the figure so the artifact says what produced it.
const N      = 40
const N_MCS  = 100
const SEEDS  = (42, 43, 44)
const ORDER  = collect(0:7)          # identity; the permuted orderings are the
                                     # test's business, not the figure's
const T_CPM  = 5.0f0
const DH_EVERY = 10                  # sweeps sampled for the ΔH pool

function measure(seed)
    rates = Float64[]
    dhs   = Float32[]
    cb = (mcs, st, dh) -> begin
        ev = count(!=(0), st)
        push!(rates, ev == 0 ? NaN : count(==(UInt8(2)), st) / ev)
        mcs % DH_EVERY == 1 && append!(dhs, dh[st .!= 0])
    end
    rp = RP(; Nr = 40, Ddot_R = 1.0, c_ext = 1.0, basis_gate_ack = true)
    RUN(; seed, n_mcs = N_MCS, N, rp, verbose = false,
        color_order = ORDER, on_sweep = cb)
    return rates, dhs
end

results = [(s, measure(s)...) for s in SEEDS]

set_theme!(Theme(fontsize = 13,
                 Axis = (spinewidth = 1.2, xgridvisible = false,
                         ygridvisible = false, xtickalign = 1, ytickalign = 1),
                 Lines = (linewidth = 2.0,),
                 Legend = (framevisible = false,)))

fig = Figure(size = (980, 400), figure_padding = (12, 20, 10, 10))

ax1 = Axis(fig[1, 1], xlabel = "Monte Carlo Steps",
           ylabel = "accepted / evaluated proposals",
           title = "Metropolis acceptance rate per sweep")
for (seed, rates, _) in results
    lines!(ax1, 1:length(rates), rates, label = "seed $seed")
end
axislegend(ax1, position = :rt)

allrates = vcat([r for (_, r, _) in results]...)
allrates = filter(!isnan, allrates)
alldh    = vcat([d for (_, _, d) in results]...)

ax2 = Axis(fig[1, 2], xlabel = "ΔH  (Metropolis acceptance functional, a.u.)",
           ylabel = "evaluated proposals",
           yscale = log10,
           title = "ΔH over evaluated proposals")
lo, hi = quantile(alldh, 0.001), quantile(alldh, 0.999)
hist!(ax2, clamp.(alldh, lo, hi); bins = 120, color = (:steelblue, 0.75))
vlines!(ax2, [0.0]; color = :black, linestyle = :dash, linewidth = 1.5)
vlines!(ax2, [Float64(T_CPM)]; color = :firebrick, linestyle = :dot, linewidth = 2)
text!(ax2, Float64(T_CPM), 1.0; text = "  T_cpm = $(T_CPM)", color = :firebrick,
      fontsize = 11, align = (:left, :bottom), space = :data)

# The provenance line is the point of the .txt sidecar: the test asserts over
# three seeds AND three colour orderings, so a figure that does not say which
# run it shows is about a different thing than the test.
Label(fig[2, 1:2],
      @sprintf("JACC port, threads backend, %d thread(s) | N=%d, %d MCS, seeds %s, color_order=identity | pooled rate %.4f (min %.4f, max %.4f), n_ΔH=%d sampled every %d sweeps",
               Threads.nthreads(), N, N_MCS, join(SEEDS, ","),
               mean(allrates), minimum(allrates), maximum(allrates),
               length(alldh), DH_EVERY);
      fontsize = 10, color = :gray30, halign = :left, tellwidth = false)

mkpath(OUT)
base = joinpath(OUT, "fig5_acceptance_regime")
save(base * ".pdf", fig)
save(base * ".png", fig; px_per_unit = 2)
@printf("wrote %s.{pdf,png}\n  pooled rate %.5f  min %.5f  max %.5f  n_ΔH %d\n",
        base, mean(allrates), minimum(allrates), maximum(allrates), length(alldh))
