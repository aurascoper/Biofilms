# Does `vols` still equal the lattice's own id histogram after a colour pass on device?
#
# This is the decisive test for lost atomic increments. `cpm_color!` maintains `vols`
# incrementally through `JACC.@atomic vols[σ] -= 1 / += 1` while writing the new label
# into `lat`. The two are redundant: a positive-id histogram of `lat` must reproduce
# `vols` exactly, on every backend, after every synchronized pass.
#
# If increments are dropped on device, volumes drift toward their starting values, the
# volume penalty stays small, and acceptance stays high -- which would be a large effect
# rather than the 5% that maximal read staleness accounts for.
#
# Launch parameters are copied verbatim from run_coupled's production call site
# (biofilms_potts_jacc.jl:536-540) so this measures the shipped kernel, not a variant.
#
# usage: julia --project=. diagnostics/volume_conflict/histogram_check.jl [N] [n_mcs]

using Printf
include(joinpath(@__DIR__, "..", "..", "biofilms_potts_jacc.jl"))

function id_histogram_host(lat, nlabels)
    h = zeros(Int32, nlabels)
    for v in lat
        v > 0 && v <= nlabels && (h[v] += Int32(1))
    end
    return h
end

function check(; N::Int = 20, n_cells_per_species::Int = 2, seed::Int = 42, n_mcs::Int = 5,
               T_cpm = 5.0f0, λ_V = 10.0f0, V_target = Int32(120))
    lat_h, spec_h, vols_h, rad_h, mel_h, _ = init_host(N, n_cells_per_species, seed, 1.0, 2.0, 1.0f0)

    lat  = JACC.array(lat_h);  vols = JACC.array(vols_h); spec = JACC.array(spec_h)
    J    = JACC.array(build_J_matrix())
    βv   = JACC.array(BETA_ION); melc = JACC.array(MEL_COEF)
    rad  = JACC.array(rad_h);  mel  = JACC.array(mel_h)
    st   = JACC.array(zeros(UInt8, N, N, N)); dh = JACC.array(zeros(Float32, N, N, N))
    Nh = N ÷ 2
    gseed = splitmix64(UInt64(seed))

    @printf("backend = %s   N = %d   parcels = %d   n_mcs = %d\n\n",
            string(JACC.backend), N, length(vols_h), n_mcs)
    @printf("%5s %7s  %10s %10s %10s  %s\n",
            "mcs", "colour", "sum(vols)", "sum(hist)", "maxdiff", "vols == histogram")

    bad = 0
    for mcs in 1:n_mcs
        fill!(st, UInt8(0)); fill!(dh, 0.0f0)
        for c in 0:7
            JACC.parallel_for((Nh, Nh, Nh), cpm_color!,
                lat, vols, spec, J, βv, melc, rad, mel,
                Int32(c & 1), Int32((c >> 1) & 1), Int32((c >> 2) & 1),
                gseed, UInt64(mcs * 8 + c), Int32(N), λ_V, V_target, T_cpm,
                st, dh)
            # to_host synchronises; the comparison is of settled state.
            lh = JACC.to_host(lat); vh = JACC.to_host(vols)
            hist = id_histogram_host(lh, length(vh))
            diffs = vh .- hist
            md = maximum(abs.(diffs))
            ok = md == 0
            ok || (bad += 1)
            if !ok || (mcs == 1) || (mcs == n_mcs && c == 7)
                @printf("%5d %7d  %10d %10d %10d  %s\n",
                        mcs, c, sum(vh), sum(hist), md, ok ? "yes" : "NO")
                if !ok
                    worst = argmax(abs.(diffs))
                    @printf("        parcel %d: vols=%d histogram=%d  drift=%+d\n",
                            worst, vh[worst], hist[worst], diffs[worst])
                end
            end
        end
    end
    println()
    @printf("passes checked = %d, mismatched = %d\n", n_mcs * 8, bad)
    return bad
end

N     = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 20
n_mcs = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 5
exit(check(; N = N, n_mcs = n_mcs) == 0 ? 0 : 1)
