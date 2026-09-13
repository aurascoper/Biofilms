# Dump the kernel's own instrumentation for ONE colour pass from an identical initial
# state, so two backends can be compared before any divergence accumulates.
#
# st: 0 never proposed, 1 evaluated-rejected, 2 evaluated-accepted
# dh: ΔH at each proposed site
#
# Comparing these separates three candidate causes that a rate number cannot:
#   proposal sets differ  -> the counter-based RNG or nb26 differ across backends
#   ΔH differs            -> delta_H arithmetic differs (Float32 ops, fma, ordering)
#   only decisions differ -> exp() or u01() differ at the accept branch
#
# usage: julia --project=. first_pass_dump.jl <out.bin>
using Printf
include(joinpath(@__DIR__, "..", "..", "biofilms_potts_jacc.jl"))

function dump_first_pass(path; N = 20, n_cells_per_species = 2, seed = 42,
                         T_cpm = 5.0f0, λ_V = 10.0f0, V_target = Int32(120), colour = 0)
    lat_h, spec_h, vols_h, rad_h, mel_h, _ = init_host(N, n_cells_per_species, seed, 1.0, 2.0, 1.0f0)
    lat  = JACC.array(lat_h);  vols = JACC.array(vols_h); spec = JACC.array(spec_h)
    J    = JACC.array(build_J_matrix())
    βv   = JACC.array(BETA_ION); melc = JACC.array(MEL_COEF)
    rad  = JACC.array(rad_h);  mel  = JACC.array(mel_h)
    st   = JACC.array(zeros(UInt8, N, N, N)); dh = JACC.array(zeros(Float32, N, N, N))
    Nh = N ÷ 2
    gseed = splitmix64(UInt64(seed))

    JACC.parallel_for((Nh, Nh, Nh), cpm_color!,
        lat, vols, spec, J, βv, melc, rad, mel,
        Int32(colour & 1), Int32((colour >> 1) & 1), Int32((colour >> 2) & 1),
        gseed, UInt64(1 * 8 + colour), Int32(N), λ_V, V_target, T_cpm, st, dh)

    sth = JACC.to_host(st); dhh = JACC.to_host(dh); vh = JACC.to_host(vols)
    open(path, "w") do io
        write(io, Int32(N)); write(io, Int32(length(vh)))
        write(io, sth); write(io, dhh); write(io, vh)
    end
    @printf("backend=%s  proposed=%d  accepted=%d  rejected=%d  sum(vols)=%d -> %s\n",
            string(JACC.backend), count(!=(0), sth), count(==(UInt8(2)), sth),
            count(==(UInt8(1)), sth), sum(vh), path)
end

dump_first_pass(ARGS[1])
