# Producer for the field-kernel timing table in docs/metal_feasibility_macos_arm64.md.
#
# Times exactly what that table describes and nothing more: the two parallel_for
# field kernels (melanin_k!, nutrient_k!) at N = 40, 200 iterations after a warm-up,
# with to_host before the clock stops so the device is synchronised. It does NOT
# launch cpm_color!, run_coupled or the host radiolysis step, so its number is a
# two-kernel microbenchmark, not the port's per-MCS cost. The backend is whatever
# LocalPreferences.toml names (JACC.set_backend); the thread count is JULIA_NUM_THREADS.
#
#   julia --project=. diagnostics/metal_feasibility/bench_field_kernels.jl [iterations]
#
# Prints one line: backend, nthreads, N, iterations, ms per step (both kernels).
include(joinpath(@__DIR__, "..", "..", "biofilms_potts_jacc.jl"))

const N_BENCH = 40
const ITERS = isempty(ARGS) ? 200 : parse(Int, ARGS[1])
const WARMUP = 20

function bench(N::Int, iters::Int)
    lat_h, spec_h, _, rad_h, mel_h, nut_h = init_host(N, 6, 42, 1.0, 2.0, 1.0f0)
    for idx in eachindex(mel_h); mel_h[idx] = Float32(0.01 * (idx % 7)); end
    mel = JACC.array(mel_h); mel2 = JACC.array(zeros(Float32, N, N, N))
    nut = JACC.array(nut_h); nut2 = JACC.array(zeros(Float32, N, N, N))
    lat = JACC.array(lat_h); spec = JACC.array(spec_h)
    αv = JACC.array(ALPHA_M); rad = JACC.array(rad_h); upt = JACC.array(UPTAKE)
    step!() = begin
        JACC.parallel_for((N, N, N), melanin_k!, mel, mel2, lat, spec, αv, rad, Int32(N), 0.5f0, 0.1f0)
        JACC.parallel_for((N, N, N), nutrient_k!, nut, nut2, lat, spec, upt, Int32(N), 0.5f0, 0.2f0, Float32(0.9))
    end
    for _ in 1:WARMUP; step!(); end
    JACC.to_host(nut2)
    t0 = time_ns()
    for _ in 1:iters; step!(); end
    JACC.to_host(mel2); JACC.to_host(nut2)   # synchronise before stopping the clock
    return (time_ns() - t0) / 1e6 / iters
end

ms = bench(N_BENCH, ITERS)
@printf("backend=%s nthreads=%d N=%d iterations=%d per_step_ms=%.3f\n",
        string(JACC.backend), Threads.nthreads(), N_BENCH, ITERS, ms)
