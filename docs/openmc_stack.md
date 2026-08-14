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

Archive SHA-256 and the `cross_sections.xml` checksum are recorded in
`~/openmc-data/*.sha256` and echoed into every `transport_result.h5`
(`nuclear_data` attr). Data-gated tests skip — reported as SKIPPED, never
passed — when `OPENMC_CROSS_SECTIONS` is unset or the file is missing.

## Test tiers

| tier | needs | gate |
|---|---|---|
| unit | numpy, h5py, pytest (bare venv) | none |
| `openmc_api` | OpenMC Python API importable | `importorskip` |
| `openmc` | executable + `OPENMC_CROSS_SECTIONS` | env-var skipif |
| Julia interop | `julia` on PATH | skipif |

Accounting is always passed / failed / SKIPPED — no `-m` deselection, so
nothing hides from the report.
