# Reference D — institutional biosafety submission, draft scope

**Status:** DRAFT · **Not an approval** · **Authorises nothing** ·
**Working table:** `data/calibration/reference_d_condition_proposal.csv`
(`authoritative_for_campaign = false`)

This is preparatory material for an institutional biosafety submission. It is not
evidence, it is not an approval, and nothing in it can clear `D-APPROVAL`. The
gate reads `data/calibration/baseline_condition.csv` by name and never this
document or its working table.

## Why the packet must freeze more than a species list

Biosafety follows **strains, not species**. Several components of the target
consortium are commonly BSL-1 while *C. neoformans* strains are commonly BSL-2,
and collection material under one species name can carry different assigned
levels. Current guidance is protocol-driven risk assessment rather than a level
label, so the chain is:

> exact strain IDs → institutional review → containment decision → approved
> procedures

An approval bound only to a species list would not actually cover the procedures
performed, which is the failure the scope digest in `approval.py` exists to make
visible after the fact.

## Scope the submission must bind

**Organisms**
- exact strain designations and their sources (collection, accession, isolate)
- containment level per strain *and per procedure*, not one level for the consortium

**Culture**
- growth medium; temperature; pH; oxygen condition
- substrate or membrane; biofilm age; flow or static operation
- culture volumes and vessels
- irradiated versus non-irradiated conditions

**Facility and apparatus**
- containment facility and suite
- irradiation source interface, if any, and how organisms and source are separated
- aerosol-generating procedures

**Handling**
- storage and transport
- decontamination and waste handling
- spill and exposure response
- training and occupational-health requirements

**Documentation**
- risk-assessment reference the approval rests on
- protocol version the approval is issued against

## What the approved row must carry back

When the approval is issued, copy the **exact approved condition** — not a
paraphrase — into `baseline_condition.csv`, with:

| field | note |
|---|---|
| `institutional_approval_id` | the identifier as issued |
| `institutional_approval_authority` | the committee or office that issued it |
| `approval_source_id` | resolves in `data/calibration/spatial/sources.csv`, `document_type = institutional_biosafety_approval`, with a url or sha256 |
| `approval_effective_date` | ISO; must not be after culturing starts |
| `approval_expiration_date` | ISO; must be after culturing starts and after today |
| `approved_protocol_version` | the version reviewed |
| `culturing_start_date` | ISO |
| `approval_scope_hash` | `approval.scope_digest(row)` over the approved conditions |
| `is_target_system` | `true` — a surrogate can never clear this gate |

`approval_scope_hash` is what makes a later edit visible: change a growth
condition without recomputing it and the gate refuses, because the approval no
longer describes the row. It detects **drift, not forgery** — recomputing it
after an edit defeats it, which is why the external document remains the anchor.

## Exit criteria

The named milestone is `REFERENCE_D_CAMPAIGN_INSTITUTIONALLY_AUTHORIZED`, and
`calibration/scripts/reference_d_status.py` prints each criterion with its state
and the date it judged against:

1. exact strain identities frozen
2. exact growth conditions frozen
3. containment facility identified
4. risk assessment referenced
5. institutional biosafety approval issued, with an issuing authority
6. approval scope matches the strains, procedures, facility and conditions
7. approval effective before culturing began, and not expired
8. approval provenance recorded and retrievable
9. `baseline_condition.csv` carries the approved target row

`CAMPAIGN_READY` requires this milestone, so the two cannot disagree. A tenth
criterion of the form "the status script agrees" would be circular and is
deliberately absent.

## What clearing it would and would not mean

It would mean exactly one thing: **the institution has authorised the specified
culturing campaign to begin.**

It would not mean Reference D has been calibrated, that a transport
configuration exists, that a pitch has been selected, or that either feedback
gate has passed. Those remain 19 unsatisfied requirements, two `PROVISIONAL`
gates, an unpopulated voxel binding, and an unset δ.

## What may proceed now

Everything that does not culture the target: synthetic fixtures, public
surrogates with `is_target_system = false`, the transport and estimator work,
the viewer stack, the rasterization ladders, and this submission packet. A
surrogate campaign is a **separate reference system** and can never clear
`D-APPROVAL` for the seven-species consortium.

δ is independent of all of this — a prospective modelling and decision-policy
declaration, not an institutional authorisation, and it must not be selected
from any effect this repository has measured.
