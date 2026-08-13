# ADDER applicability decision

**Decision: ADDER is not in this project's runtime path.** This document is
the decision gate the review required in place of treating ADDER as a
drop-in OpenMC material-update layer.

## What stock ADDER is

Argonne's ADDER (v1.1.0 user guide) is a reactor fuel-management and
depletion system. Its supported external codes are **MCNP5/6** for transport
(tally extraction via MCNPTools) and **CRAM / ORIGEN2** for depletion and
decay. It has an abstract neutronics interface designed for extension, so an
OpenMC backend is architecturally conceivable — but that means implementing
model parsing, material synchronization, execution, flux/reaction-rate
extraction, volume handling, solver results, and the supported
fuel-management operations. That is a separate research-software project,
not a thin wrapper around a mesh tally.

## Why it does not fit this repository's physics

1. **The repository's radiation field is gamma-only.** OpenMC's photon
   transport does not implement photoneutron reactions, so a gamma-only
   source establishes no transport mechanism for neutron-driven thorium
   depletion. The draft's thorium half-life discussion does not supply a
   spectrum, activity, or decay-chain equilibrium assumption — those are
   REQUIRED config inputs, not derivable.
2. **The updates the draft assigned to ADDER are not depletion.** Biological
   aging, polymer/membrane damage, biosorption, and cell movement are
   material-chemistry and biological-state updates. They belong to the
   CPM/radiodialysis layer (and the config-defined material mapper), not to
   a nuclide-inventory solver.

## Decision matrix

| Physical case | Recommended runtime |
|---|---|
| External gamma field; local dose + polymer/biological response (THIS REPO) | OpenMC photon transport + the existing/improved chemistry-damage model. **No depletion solver.** |
| Neutron field with activation/transmutation | OpenMC transport + `openmc.deplete` |
| ADDER specifically required (reactor fuel management) | MCNP + ADDER |
| Research goal: add OpenMC support to ADDER | Separate ADDER fork/backend project |

## Revisit triggers

Reopen this decision only if (a) a neutron source term enters the problem
definition with a declared spectrum and rate, or (b) isotope inventory
evolution (activation, transmutation, decay-chain buildup) becomes a
scientific question — in which case the first stop is `openmc.deplete`,
not ADDER.
