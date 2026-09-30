# Ellen review package, 2026-09-30

This is the repository-built review draft. It is separate from the CopyDesk PDF and from the stable 2026-09-29 review artifact.

## Included work

- completed coefficient, 400-MCS horizon and same-species-contact follow-ups;
- corrected synthetic-survival specification with separate survival, survivor-radius and identified-parcel-displacement outcomes;
- a Dsup-inspired protection contrast at fixed absorbed dose, including the ortholog-specific bacterial result and the reason it does not supply a generic protection coefficient, with repair and dormancy deferred to their own state models;
- a read-only phase passenger specification and its open evidence gates;
- a four-arm design for testing synthetic survival and the phase passenger alone
  and together, with the first combined model limited to removal of nonviable
  parcels from the active phase graph;
- a public-data ART and microbially influenced corrosion study plan; and
- ART motif and mechanism questions for Ellen.

The survival, protection, phase and four-arm sections are prospective
specifications. No hazard, protection, repair, dormancy, phase or factorial
result is reported.

## Validation

- the full Julia test suite exited successfully;
- all 675 Python calibration, manuscript, claims, figure, vocabulary and scope guards passed;
- three consecutive `pdflatex` passes completed without undefined references; and
- Figure 1 passed its render-staleness guards, and the revised result, future-work, survival, phase and reference pages were inspected after rendering.

The repository remains private. Access is available on request.
