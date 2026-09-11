# Metal feasibility on Apple Silicon — three gates, reported including failures

First evaluation of the JACC port on a Metal backend. Apple M4 (4 performance + 6
efficiency cores), macOS 26.6.2, Julia 1.12.6 aarch64, JACC 1.3.1, Metal.jl 1.11.0,
device `AGXG16GDevice`. The companion result for AMD is
[`jacc_coupling_port.md`](jacc_coupling_port.md); this file is deliberately its sibling
and reaches a compatible conclusion by a different route.

**Verdict: do not port.** Metal runs the port correctly and is 2.3x slower than the
CPU threads backend at saturation. That confirms what the ROCm measurement already
said — the near-term work is kernel structure, not more FLOPs.

`Metal` was added to `Project.toml` and `AMDGPU` removed for the duration of these
runs; both are local and uncommitted, and `Project.toml` is byte-identical to HEAD in
this commit. `JACC.set_backend` writes `LocalPreferences.toml`, which `.gitignore:379`
already excludes.

## Gate 1 — precision: FAILS, and blocks nothing

```
Metal.functional() : true      device : AGXG16GDevice
Float32            : OK        (MtlVector{Float32, PrivateStorage})
Float64            : ErrorException("Metal does not support Float64 values, try using Float32 instead")
Int64              : OK        (MtlVector{Int64, PrivateStorage})
```

The refusal reaches through JACC unchanged: `JACC.array(ones(Float64, 8))` throws the
same error with `array_type() == Metal.MtlArray`. **JACC does not paper over a hardware
limit**, and it should not be expected to.

**This blocks nothing, and the reason corrects a standing assumption.** The framing
this gate was written against was: the CPM is `Float64` throughout, so a Metal port
*would become* a `Float32` model and could not be bit-identical. The conclusion is
right and the tense is wrong. The **serial** CPM is `Float64`. The **JACC port never
was** — `biofilms_potts_jacc.jl` is `Float32` end to end (`build_J_matrix` returns
`J32`, the kernels take `0.5f0`/`0.1f0`/`Float32(0.9)`, the selftest allocates
`zeros(Float32, N, N, N)` and compares against serial at `rtol = atol = 1e-4`).

That decision was made when the port was written, not deferred to this gate. It is
also already reflected in the suite: `tests/jacc_port_tests.jl` carries `Float32`
tolerances, and `tests/contract_csv.jl` does not and must never be run against a GPU
backend.

## Gate 2 — backend: PASS

`JACC.set_backend("metal")`, restart, `biofilms_potts_jacc.jl --selftest`:

```
[ Info: Metal backend loaded
selftest: delta_H matches serial on 200 site pairs
selftest: melanin/nutrient kernels match serial one-step update
selftest: PASS
```

Verified to be a real device path rather than a silent fallback:
`JACC.backend == metal`, `JACC.array_type() == Metal.MtlArray`, and
`JACC.array(ones(Float32, 8))` returns `MtlVector{Float32, PrivateStorage}`. The repo
uses only `array`, `array_type`, `backend`, `parallel_for`, `to_host` and
`set_backend`, all of which `JACC 1.3.1`'s `ext/MetalExt` provides.

**One probing trap worth recording.** `using JACC` alone leaves the backend
unregistered — `JACC.array_type()` then throws `MethodError: no method matching
get_backend(::Val{:metal})`, which reads exactly like "Metal is broken." It is not:
`MetalExt` is a package extension, and `JACC.@init_backend` (what
`biofilms_potts_jacc.jl:25` calls) is what loads it. A probe that omits the macro
measures its own omission.

## Gate 3 — agreement: as specified it PASSES BUT CANNOT FAIL

The gate was specified as `validate_serial.jl 42`, with the note that it may
legitimately fail on `Float32`. It passes — **8/8 CSV rows byte-identical** to
`tests/fixtures/serial_seed42.csv` with the Metal preference active.

**It could not have done anything else.** `validate_serial.jl` contains no reference
to `JACC`, `Metal` or any backend; it exercises the `Float64` serial monolith. The
backend preference cannot reach it, so it can detect nothing whatever about the Metal
port. The pass is real and it is evidence for exactly one proposition: selecting a GPU
backend does not disturb the serial path.

The agreement test that can fail is `tests/jacc_port_tests.jl`, which is the one
carrying `Float32` tolerances:

```
[ Info: Metal backend loaded
[ Info: JACC port selftest running  backend = "MetalBackend"
Test Summary:                                 | Pass  Total   Time
JACC port kernels versus the serial reference |    9      9  13.9s
```

**Real Gate 3: PASS, 9/9 under `MetalBackend`.**

## Timing, and the comparison that is honest

Both `parallel_for` kernels (`melanin_k!`, `nutrient_k!`) at N = 40, 200 iterations,
warmed up, `to_host` before stopping the clock so the device is synchronised:

| Backend | per step | vs saturated CPU |
|---|---|---|
| `metal` | **0.739 ms** | **2.3x slower** |
| `threads`, nthreads = 1 | 20.199 ms | 63x slower |
| `threads`, nthreads = 2 | 0.623 ms | 1.96x slower |
| `threads`, nthreads = 4 | **0.318 ms** | saturated (4 P-cores) |
| `threads`, nthreads = 8 | 0.326 ms | no further gain |
| `threads`, nthreads = 10 | 0.318 ms | no further gain |

**Metal beats nthreads = 1 by 27x, and that comparison is worthless.** The 1-thread
figure is not this kernel's serial speed. `JACC/src/threads/threads.jl:12-20` branches
on `Threads.nthreads() == 1` and runs the bare loop, taking `Polyester.@batch` only
above one thread. The 1→2 step is therefore 32x — Polyester versus a plain closure
loop, not parallelism. Real parallel scaling is the 2→4 step, a clean 2x that
saturates at the 4 performance cores and is flat thereafter.

Against the saturated CPU — the only honest comparison — **Metal is 2.3x slower**,
and the reading is the same one `jacc_coupling_port.md` gives for ROCm at 8x: at
N = 40 the lattice is ~6.4e4 sites and the run is launch-overhead-bound, not
arithmetic-bound. A larger device does not fix that; kernel structure does.

**How to read the nthreads = 1 row.** It is not a defect in this repository.
`tests/jacc_parity_tests.jl:235` asserts `Threads.nthreads() == 1` deliberately,
because its byte-identical table comparison is interpretable only on the single-thread
configuration it was measured on — with four threads, scheduling changes the tables
even when acceptance does not read the nutrient field. That is a correctness pin, and
it costs 63x throughput in the JACC threads path. Worth knowing before anyone reads a
JACC threads timing taken at the default and concludes the CPU path is slow.
`JULIA_NUM_THREADS` is not exported anywhere in this environment, so 1 is what any
unconfigured invocation gets.

## The ANE

The Apple Neural Engine is reachable only through CoreML. There is no Julia path, no
JACC backend, and nothing about a Metropolis update on a Potts lattice that CoreML
expresses. It is unavailable for this model. No hedge.

## Reproducing

```
julia --project=. -e 'using Pkg; Pkg.add("Metal")'          # and remove AMDGPU: no darwin artifacts
julia --project=. -e 'using JACC; JACC.set_backend("metal")'
julia --project=. biofilms_potts_jacc.jl --selftest          # Gate 2
julia --project=. validate_serial.jl 42                      # Gate 3 as specified
julia --project=. -e 'using JACC; JACC.set_backend("threads")'   # restore
```

## Observed, not fixed

`libopenvkl_module_cpu_device.dylib` is absent from the macOS ParaView bundle, so
every `pvpython` run prints `[openvkl] INITIALIZATION ERROR`. Unrelated to JACC;
recorded in `docs/visualization/macos-arm64-2026-09-11/README.md` and repeated here
only because both macOS findings landed the same day.
