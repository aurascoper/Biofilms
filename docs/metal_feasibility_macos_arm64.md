# Metal feasibility on Apple Silicon — three gates, reported including failures

First evaluation of the JACC port on a Metal backend. Apple M4 (4 performance + 6
efficiency cores), macOS 26.6.2, Julia 1.12.6 aarch64, JACC 1.3.1, Metal.jl 1.11.0,
device `AGXG16GDevice`. The companion result for AMD is
[`jacc_coupling_port.md`](jacc_coupling_port.md); this file is deliberately its sibling
and reaches a compatible conclusion by a different route.

**Verdict: do not port — and the reason is correctness, not speed.** The CPM kernel executes
on device, but the port does **not reproduce the CPU's acceptance statistics**: the Metropolis
acceptance rate is **1.94x** the threads-backend rate, uniformly across all eight parity
classes. Timing is the lesser finding.

An earlier revision of this file said "Metal runs the port correctly and is 2.3x slower."
That was measured over **two field kernels** (`melanin_k!`, `nutrient_k!`) and a host-side
`delta_H`; `selftest` never launches `cpm_color!`, so it could not and did not establish
whole-port correctness. Running the parity tier on device is what found the defect.

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
right and the tense is wrong. The **serial** CPM is `Float64`; the port's **device-facing
CPM and field arrays never were** (`build_J_matrix` returns `J32`, the kernels take
`0.5f0`/`0.1f0`/`Float32(0.9)`, the selftest allocates `zeros(Float32, N, N, N)` and compares
against serial at `rtol = atol = 1e-4`).

**"`Float32` end to end" is wrong and an earlier revision of this paragraph said it.** The
same file's host radiolysis is `Float64` — see the inventory below. The accurate statement is
that the *device-resident arrays* are `Float32`, which is why Metal's refusal blocks nothing;
it is not a claim about the file.

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

## The CPM kernel on device: 51 pass, 14 fail

`tests/jacc_parity_tests.jl` is the tier that drives `run_coupled -> cpm_color!`, including
the `JACC.@atomic` volume updates at `biofilms_potts_jacc.jl:217-218`. Run through
`tests/runtests.jl` with `backend = metal`:

| Testset | Metal |
|---|---|
| the class encode/decode pair agrees with the kernel | **16/16 pass** |
| `color_order` must be a permutation of 0:7 | 3/3 pass (host-side) |
| per-sweep reset, discriminator, all eight classes | **2/2 pass** — `nonfinite == 0` |
| the guard can fire | **2/2 pass** |
| no decomposition artifact across seeds and orderings | 24 pass, **12 fail** |
| the exemption's claim (byte-identical tables) | 4 pass, **2 fail** |

Same code, same `Project.toml`, only `LocalPreferences.toml` changed:

| | acceptance rate over 3 seeds x 3 orderings | band (0.14, 0.23) |
|---|---|---|
| `threads` | 0.17061 .. 0.19137, mean **0.17933** | all inside; 65/65 pass |
| `metal` | 0.33121 .. 0.36889, mean **0.34870** | all outside |

**1.944x.** Not a marginal band excursion.

### The failure breakdown, from the raw log

An earlier revision of this file said "every failure in that testset is the rate". **That was
wrong**, and the arithmetic was checkable without the log: nine runs carry four assertions
each, so 24 passes and 12 failures cannot be nine rate failures alone. Counting only the
`Test Failed at` lines in the run:

| assertion | line | failures |
|---|---|---|
| `RATE_BAND[1] < s.rate < RATE_BAND[2]` | 200 | **9** of 9 |
| `s.V < V_MAX` (pooled) | 201 | 0 |
| `s.maxdev < MAXDEV_MAX` (pooled) | 202 | 0 |
| `s1.V < V_MAX && s2.V < V_MAX` (half-windows) | 209 | **3** of 9 |
| testset 6 table comparisons | 250, 251 | 1 each |
| `Threads.nthreads() == 1` | 235 | **0** — passed, as predicted |

**So a time-dependent parity effect is present and was wrongly excluded.** The pooled
statistics do pass everywhere — Cramer's V 0.00873..0.02382 against `V_MAX = 0.025`, max
per-class deviation 0.0223..0.0503 against `MAXDEV_MAX = 0.12`. But in **3 of 9 runs** a
half-window V exceeds `V_MAX`: within a time window the colour classes do **not** agree, and
pooling over the full run averages that away.

Julia prints no operand values for a `&&` compound, so the specific seeds, orderings and
halves are not recoverable from this log. Re-running with the two halves asserted separately
is required before any statement about which window drifts.

The global rate inflation is real and large. It is **not** established that it is uniform
across classes, and "not a decomposition artifact" is not supported by this run.

### The mechanism, stated as a reading of the source

`cpm_color!` is order-independent in three of its four state interactions. The RNG is
counter-based on `(seed, step, linear index)`, so draws do not depend on execution order. The
checkerboard makes concurrently-updated sites spatially non-adjacent, so `lat` reads are
race-free — the kernel's own comment says so. `dh` and `st` are per-site write-only.

**`vols` is the exception.** `delta_H` *reads* `vols` for the volume-constraint term while
`JACC.@atomic vols[sigma] -= 1` / `+= 1` mutate it. Two sites in one colour class are
guaranteed spatially non-adjacent but **not** guaranteed to belong to different cell labels, so
whether one site's ΔH sees another's volume update depends on interleaving. Atomic addition
commutes for the final total; it does not make the intermediate reads deterministic. **The direction of the bias is not determined by undercounting alone.** For a copy from
positive label s into positive label t the volume contribution is

    ΔH_V = λ_V[2(V_s − V_target) + 1] + λ_V[−2(V_t − V_target) + 1]

so with read errors e_s and e_t relative to a sequential reference the energy error is
**2λ_V(e_s − e_t)**. Underestimating the **donor** volume lowers ΔH; underestimating the
**recipient** volume raises it. An earlier revision of this file asserted that undercounting
lowers ΔH and therefore accepts more often. That is false as stated — the two terms carry
opposite signs, and the net sign must be **measured**, not inferred.

This is a reading of the source consistent with the measurement. It is **not** a proof that the
`vols` read is the only contributor, nor that it accounts for the factor 1.944, and no attempt
was made here to isolate it.

### Testset 6 was declared uninterpretable before the run, and is

`jacc_parity_tests.jl:235` asserts `Threads.nthreads() == 1` because its byte-identical table
comparison "is interpretable only on the single-thread configuration it was measured on". That
assertion **passes on Metal** — there is one Julia host thread — while thousands of work-items
execute concurrently. One host thread does not serialize GPU execution.

Its 2 failures are tables differing in single counts (76 vs 75, 117 vs 116). That is
non-determinism, and it says nothing either way about the testset's actual claim — that
acceptance does not read the gated biomass basis. **The result is recorded and is not evidence
about the exemption.** The claim remains established on the CPU tier, which is unchanged.

### What this does not establish

No full-workload timing was produced. The CPU tests are unrelaxed and still pass 65/65. This
says nothing about ROCm, where the decomposition was separately measured.

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

## The Float64 inventory is wider than previously recorded

`docs/jacc_coupling_port.md:133-139` concedes a Float64 parameter side but names only
`CPMParams`. Two more live in **the port's own file**: `RadiolysisParams`
(`biofilms_potts_jacc.jl:266-288`, 12 of 14 fields `Float64`) and `RadiolysisState`
(`:290-304`, every numeric field `Float64`, including three `Vector{Float64}`), stepped on host
every MCS from `run_coupled`. "The JACC port is Float32 throughout" is wrong; the device-facing
arrays are Float32 and the host radiolysis integrator is not. Metal's `Float64` refusal does
not force their conversion unless they move to the device.

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
