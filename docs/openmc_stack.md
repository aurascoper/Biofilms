# Pinned OpenMC stack

Per review amendment 10: no floating `pip install openmc`. The transport stack
is a pinned conda-forge environment; every result file records the versions,
seed/particles/batches, and the nuclear-data identity used.

## Environment

Created with micromamba (conda-forge):

```
micromamba create -n openmc-biofilms -c conda-forge "openmc=0.15.3" h5py numpy pytest
```

Resolved on the development machine (2026-08-13):

| package | version | build |
|---|---|---|
| openmc | 0.15.3 | dagmc_nompi_py313h6424856_102 |
| python | 3.13.15 | hb101c97_101_cp313 |
| hdf5 | 1.14.6 | nompi_h19486de_110 |
| h5py | 3.16.0 | nompi_py313h253c126_102 |
| numpy | 2.5.2 | py313hf6604e3_0 |

The conda package ships both the `openmc` executable and the Python API from
the same build — no source-checkout skew.

## Nuclear data

Official OpenMC ENDF/B-VIII.0 HDF5 library (includes incident-neutron,
**photoatomic**, and atomic-relaxation data) from https://openmc.org/data/:

```
curl -L -o endfb-viii.0-hdf5.tar.xz \
  https://anl.box.com/shared/static/uhbxlrx7hvxqw27psymfbhi7bx7s6u6a.xz
tar -xJf endfb-viii.0-hdf5.tar.xz
export OPENMC_CROSS_SECTIONS=$HOME/openmc-data/endfb-viii.0-hdf5/cross_sections.xml
```

The **archive** SHA-256 is recorded in `~/openmc-data/*.sha256`. There is no
separate checksum for `cross_sections.xml`, and an earlier version of this
document claimed one — the archive digest is the stronger identity in any case,
since it pins every data file rather than only the index.

It is echoed into every `transport_result.h5` (`nuclear_data` attr). Note that
`synthetic_e2e.py` writes the **literal string** `"endfb-viii.0"` rather than
reading the environment, so it records an assumption about the library rather
than an observation of it; `openmc_nested_pilot.py:_nuclear_data_id()` resolves
`OPENMC_CROSS_SECTIONS` and its recorded digest instead, and names the file the
digest came from. Data-gated tests skip — reported as SKIPPED, never passed —
when `OPENMC_CROSS_SECTIONS` is unset or the file is missing.

## Running anything in the transport tier

The env is not on `PATH` and `OPENMC_CROSS_SECTIONS` does not persist between
shells, so every transport command is these two lines. A session that skips
them sees `NOT EVALUATED` gate verdicts and SKIPPED tests, which is the
designed behaviour and is easy to misread as "OpenMC is not installed":

```sh
export OPENMC_CROSS_SECTIONS=$HOME/openmc-data/endfb-viii.0-hdf5/cross_sections.xml
micromamba run -n openmc-biofilms <command>
```

The packages must be installed into that env too, `contract/` first because it
is a path dependency with no PyPI presence:

```sh
micromamba run -n openmc-biofilms pip install -e contract -e "coupling[dev]"
```

## Test tiers

| tier | needs | gate |
|---|---|---|
| unit | numpy, h5py, pytest (bare venv) | none |
| `openmc_api` | OpenMC Python API importable | `importorskip` |
| `openmc` | executable + `OPENMC_CROSS_SECTIONS` | env-var skipif |
| Julia interop | `julia` on PATH | skipif |

Accounting is always passed / failed / SKIPPED — no `-m` deselection, so
nothing hides from the report.

## Golden-tally fixture

`coupling/tests/fixtures/golden_tally_water_phantom.json` pins 12 REAL OpenMC
runs (2 outer draws x 3 replicates x {baseline, feedback}, water-phantom
geometry, feedback density x1.35 — the same `DENSITY_SCALE` lever
`openmc_nested_pilot.py` already uses) — the raw per-source heating tally,
not anything derived from it. `coupling/tests/test_gate_composition.py`
replays the pin through the real `specific_energy_per_source ->
debiased_squared_effect -> decide()` chain in the ordinary (no-OpenMC) test
tier, closing the one seam nothing else in this repo tests: a gate decision
made from something a tally actually produced, not from hand-fabricated
`numpy` arrays.

Regenerate with `coupling/scripts/regenerate_golden_tally.py` under this
env (needs `OPENMC_CROSS_SECTIONS`, same as everything else here). It
refuses to overwrite the committed fixture unless all 12 runs complete,
matching `openmc_nested_pilot.py`'s `writes_canonical_tables` /
`refuse_partial_publish` guard. Legitimately changes only if OpenMC or the
nuclear-data library changes — the fixture header records both — which is
what `.github/workflows/golden-tally-verification.yml` regenerates and
diffs against on exactly that trigger, same compare-only-in-CI discipline
as `tests/contract_csv.jl` for the serial fixture.
