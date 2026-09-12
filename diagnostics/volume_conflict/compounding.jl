# Does a per-pass divergence of a few percent compound to the observed 1.944x?
#
# The parity tier runs n_mcs = 50, i.e. 400 colour passes. A single pass from a common
# state is NOT the right comparison for a 50-MCS steady-state acceptance rate, and an
# earlier revision of this diagnostic made exactly that error: it measured one pass,
# found ~5%, and concluded the volume conflict could not explain a 94% gap.
#
# Both arms here run the full trajectory. They differ only in when `vols` is read.
using Printf
include(joinpath(@__DIR__, "oracle.jl"))

function trajectory(; N = 20, n_cells = 2, seed = 42, n_mcs = 50,
                    λV = 10.0f0, Vt = Int32(120), T = 5.0f0, snapshot::Bool)
    lat, spec, vols, rad, mel, _ = init_host(N, n_cells, seed, 1.0, 2.0, 1.0f0)
    for idx in eachindex(mel); mel[idx] = Float32(0.01 * (idx % 7)); end
    J = build_J_matrix()
    proposed = 0; accepted = 0
    trace = Tuple{Int,Float64}[]
    for mcs in 1:n_mcs
        for c in 0:7
            ox, oy, oz = c & 1, (c >> 1) & 1, (c >> 2) & 1
            props = colour_pass!(lat, vols, spec, J, BETA_ION, MEL_COEF, rad, mel,
                                 ox, oy, oz, splitmix64(UInt64(seed)), UInt64(mcs * 8 + c),
                                 N, λV, Vt, T; snapshot = snapshot)
            proposed += length(props)
            accepted += count(p -> p.accepted, props)
        end
        mcs % 10 == 0 && push!(trace, (mcs, accepted / proposed))
    end
    return proposed, accepted, trace
end

function main()
    println("Both arms, full 50-MCS trajectory. Only the vols read timing differs.\n")
    @printf("%-28s %10s %10s %10s\n", "arm", "proposed", "accepted", "rate")
    pl, al, tl = trajectory(snapshot = false)
    ps, as, ts = trajectory(snapshot = true)
    @printf("%-28s %10d %10d %10.5f\n", "live vols (reference)", pl, al, al / pl)
    @printf("%-28s %10d %10d %10.5f\n", "sweep-start snapshot", ps, as, as / ps)
    @printf("%-28s %32.4fx\n", "ratio", (as / ps) / (al / pl))
    println()
    println("  cumulative rate every 10 MCS:")
    @printf("  %6s %12s %12s %10s\n", "mcs", "live", "snapshot", "ratio")
    for i in eachindex(tl)
        m, rl = tl[i]; _, rs = ts[i]
        @printf("  %6d %12.5f %12.5f %10.4f\n", m, rl, rs, rs / rl)
    end
    println()
    println("  measured on Metal vs threads at the same n_mcs: 0.34870 / 0.17933 = 1.944x")
end

main()
