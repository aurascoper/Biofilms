# ADR: an inert generic signal over the verified label trajectory

Status: implemented as a separate diagnostic; acceptance coupling is excluded.
Claim key: QS-INERT-01. Basis: declared. PR #25 remains on hold.

The new Julia module integrates `dA/dt = D_A laplacian(A) + p_s I_s - gamma A`
on the fixed interior mask, with zero signal at masked walls and no flux across
outer cube faces, as in the nutrient stencil. Production is per occupied lattice voxel. Empty interior sites diffuse
and decay signal but produce none. The initial field is exactly zero.

The existing nutrient solver supplies the structural example, but this is an
independent module with decay and positivity-preserving Euler substeps. It reads
the authoritative 101 snapshots, MCS 0–100, after verifying their pinned archive
manifest. The species snapshot at MCS k is held over [k,k+1): **left-endpoint**
occupancy. It does not infer physical arrival times within the MCS.
The signal uses its own declared MCS clock. The older `dt_field=0.5` setting is
not a measured time conversion for this new field.

No CPM state, Hamiltonian coefficient, RNG, complete restart, or gated uptake
quantity is passed to this module. The original run is never rewritten. Its
accepted-copy record remains the source of event order, and the signal remains
a sampled field. At an intermediate copy the viewer must not imply a measured
sub-MCS signal value. The serial-source lexical guard is deliberately narrower
than a general dependency proof; unchanged parent artifact hashes are the
runtime check for this diagnostic's immutable input trajectory.

## Parameter declaration

| Quantity | Demo value | Numerical unit | Basis / physical status |
|---|---:|---|---|
| D_A | 0.2 | lattice-site squared / MCS | declared / blocked |
| gamma | 0.1 | inverse MCS | declared / blocked |
| p_s, each of seven species | 1 | arbitrary signal units / occupied voxel / MCS | uniform numerical illustration / unmeasured |
| A_threshold | 5 | arbitrary signal units | declared display endpoint / no biological switch |
| initial A | 0 | arbitrary signal units | declared initial condition |
| maximum substep | 0.5 | MCS | numerical, reduced further for positivity |

Physical conversion needs D-PITCH and D-TIMESERIES and signal calibration. No
molecule is assigned to any species, no production rate is inferred from a
species name, and this common scalar does not establish communication across
the seven organisms. A molecule-by-species evidence audit is needed before
assigning biological rates, sensing specificity or receptors. Hill K_A and n
are absent: no response function exists in this diagnostic. D-APPROVAL and
delta remain unchanged and unset, respectively.

## Endpoint and geometry

The reported endpoint is `count(occupied interior AND A >= A_threshold) /
count(occupied interior)` at each saved MCS. The numerator and denominator are
both recorded. An empty denominator is missing, not zero. This is a numerical
consequence of declared transport/source assumptions, not independent evidence
of biological quorum sensing, EPS production, membranes or reconstitution.

For a **steady, uniformly producing slab** with half-width L and A(+/-L)=0,
`A(0)=p/gamma * (1-sech(L sqrt(gamma/D)))`. Solving for a threshold gives
`L=sqrt(D/gamma) acosh(1/(1-gamma A_threshold/p))`, finite only when
`0 <= gamma A_threshold < p`. At gamma=0 its limit is
`sqrt(2 D A_threshold/p)`. The tests hold the slab width fixed while refining
the grid. This half-width is not a critical radius for a transient, patchily
occupied 3D cluster. Nor does the equation guarantee that every parcel cluster
crosses the selected threshold or that its geometric centre crosses first.

## Reproduce

The first development prototype used post-copy sources and absorbing cube
faces. The delivered run uses the later supplied plan's left-endpoint sources
and Neumann cube faces; their endpoint numbers must not be interchanged. The
prototype is retained separately, not used for the delivered figure.

From the repository root, on Linux or macOS with Julia 1.12:

```sh
julia --project=diagnostics/inert_signal -e 'using Pkg; Pkg.instantiate()'
julia --project=diagnostics/inert_signal diagnostics/inert_signal/test_numerics.jl
julia --project=diagnostics/inert_signal diagnostics/inert_signal/run.jl PARENT_RUN NEW_OUTPUT EXPECTED_PARENT_MANIFEST_SHA256 diagnostics/inert_signal/demo.toml
```

NEW_OUTPUT must not exist. The derived manifest pins the input manifest, code,
Julia executable, project, parameters, field files, metrics and ParaView series.
Open `paraview/signal_trajectory.pvd`, colour cell data by species, cell_id or
signal, and animate MCS 0–100. Use categorical colours for labels and continuous
colours for signal. Ray tracing is not required. Temporal Statistics can display
binary indicator means; it cannot reconstruct ordered returns.

The archive-pinned manuscript and Figure 5 remain evidence of the original CPM
labels. New signal output is a separate diagnostic, not retroactively part of
that manuscript's model specification. A later acceptance coupling would need
another ADR, claims and parameter provenance entries, and a revised specification.

The portable browser explorer is tested with the actual data at narrow and wide
sizes in light and dark themes. Its optional accepted-copy replay was completed
before the later plan deferred overlays; it does not advance or display the
signal between saved MCS values. The native GLMakie implementation is provided
by `viewer/visualize_signal.jl`. Its HDF5 reader is tested without OpenGL; the
native window remains untested in this hosted environment and needs a target
Mac/Linux check. Run `julia --project=viewer viewer/visualize_signal.jl PARENT_RUN
DERIVED_RUN` after instantiating that viewer project.

The existing `lattice_evidence.jl verify` command rewrites verification receipts
even when the scientific checks pass. Run it on an isolated archive extraction,
not on the immutable authoritative parent used by this diagnostic.
