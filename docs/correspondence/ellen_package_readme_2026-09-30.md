# Ellen review package, 2026-09-30

This is the repository-built review draft. It is separate from the CopyDesk PDF and from the stable 2026-09-29 review artifact.

## Included work

- completed coefficient, 400-MCS horizon and same-species-contact follow-ups;
- corrected synthetic-survival specification with separate survival, survivor-radius and identified-parcel-displacement outcomes;
- a Dsup-inspired protection contrast at fixed absorbed dose, including the ortholog-specific bacterial result and the reason it does not supply a generic protection coefficient, with repair and dormancy deferred to their own state models;
- a read-only phase passenger specification and its open evidence gates;
- a four-arm design for testing synthetic survival and the phase passenger alone
  and together, with the first combined model limited to removal of nonviable
  parcels from the active phase graph, paired checks for unchanged survival and
  CPM trajectories, and analysis that separates selective removal from locking;
- a public-data ART and microbially influenced corrosion study plan; and
- ART motif and mechanism questions for Ellen.

The survival, protection, phase and four-arm sections began as prospective
specifications. The survival, passive-phase and growth parts were then
implemented and executed from source `841d052` on 2026-09-30; their results
are appended in `diagnostics/fitness_phase/PREDICTION.md`. The protection,
repair and dormancy parts remain specifications, and no result is reported
for them.

## Validation

- the preceding `2ffa9ea` build passed the Julia test entry point and all 675
  Python calibration tests;
- this specification refinement passed all 142 Julia manuscript assertions and
  135 Python manuscript, claims, citation, vocabulary and figure checks;
- three consecutive `pdflatex` passes completed without undefined references
  or overfull boxes; and
- the revised joint-exercise page and its following page were inspected after
  rendering. Figure 1 passed its render-staleness guards.

These checks cover the documents and the implementation at the time of this
package. The paired survival/phase invariants and selection controls were
later executed from `841d052`: 272 of 272 seed records are paired-exact and
restart-exact. That execution postdates this package; see
`diagnostics/fitness_phase/PREDICTION.md`.

The repository remains private. Access is available on request.
