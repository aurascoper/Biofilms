# MURR as a candidate site for the lab-gap campaign — unsourced, forward-looking

2026-08-24 · zero prior mentions in this repository (`grep -rin "MURR\|Missouri Research
Reactor\|reactor pool\|hot cell\|beamport\|Lu-177\|Lutetium\|isotope production" data/ docs/`
returns nothing) — genuinely new territory, no existing convention in this repo to match
against. This document does not create one; it records a candidate lead honestly.

**No `source_id`.** Nothing below is cited from a document in `data/sources.csv`, because
none exists yet for MURR. Nothing in this note may be treated as fact, quoted elsewhere in
this repository, or promoted into any `data/research/*_candidates.csv` row until a real
citable source — a MURR facilities/capabilities page, published documentation, or written
confirmation from a contact — exists and is pinned the way every other source in
`data/sources.csv` is (`title, publisher, document_version, url, accessed_date, sha256`).

## Why this exists

`docs/research/radiotrophic_lab_gap.md` already specifies the smallest experiment public
data cannot replace: one inoculum, one reactor run, paired sibling coupons split across
destructive assays, matched blanks, and a dosimetered irradiation — gated throughout on an
`institutional_approval_id` that only an institutional biosafety committee can issue. That
document is silent on *where* such an irradiation would happen. Hunter has a live, unrelated
introduction in progress — arranged by Caixia Wan, to a MURR-affiliated nuclear-research
contact (see `Developer/jobs/career-ops/interview-prep/murr-caixia-wan.md`, which already
frames everything below as questions to ask, not facts to present) — and MURR is a
physically plausible venue for exactly that gap. This note exists so that possibility is
written down honestly rather than either forgotten or silently treated as more solid than it
is.

## What's plausible-and-likely-true vs. unverified extrapolation

**Plausible, matches MURR's public reputation** (still no `source_id` pinned here — treat as
`unverified`, not `false`, per this repo's own boolean convention that a blank/unsourced cell
is never read as zero):
- Operates at 10 MW; described as the highest-power, highest-flux university research
  reactor in the United States.
- Operates near-continuously (close to 24/7, year-round).
- A significant producer of medical and research-grade radioisotopes (e.g. a Lu-177-style
  neutron-capture production route from a stable target), which implies real neutron-flux
  and radiochemistry capability on site.
- Physically co-located with University of Missouri biological sciences/biochemistry
  facilities.

## What's NOT licensed by anything on record

1. **Not `institutional_approval_id`.** Nothing here is or substitutes for a biosafety
   determination. Per `Biofilms/AGENTS.md` rule 6, that determination belongs to an
   institutional biosafety committee and follows strains, not species — it cannot be
   asserted by an agent, a conversation, or this document.
2. **Not a confirmed experimental service.** No source confirms MURR offers, or would grant,
   biological-sample irradiation studies of the kind `radiotrophic_lab_gap.md` specifies, to
   an outside or early-career researcher. "Irradiation chambers, beamports, and pneumatic
   tubes ... researchers could place fungal colonies into specialized capsules and lower them
   into the reactor pool or beamlines" is a plausible-sounding extrapolation from general
   knowledge of research-reactor design, not a documented MURR capability or offer.
3. **Not a dosimetry guarantee.** Even if irradiation access exists, `radiotrophic_lab_gap.md`
   requires a full dosimetered, biosafety-cleared, paired-sibling campaign design — a
   available irradiation source alone does not satisfy that; the ISS study this repo already
   audited is the standing example of an irradiation with no stated dosimeter at all.
4. **Not a target-gate clearance.** Even a real, executed MURR collaboration would still need
   to clear every gate `radiotrophic_lab_gap.md` and the compatibility audit already specify
   (paired wet/dry mass, matched blanks, closed elemental composition, a calibrated 3D
   volume) — nothing about the venue shortcuts the measurement design itself.

## STOP CONDITIONS

- **Abandon rather than cite this document as a source.** It has no `source_id` and isn't
  one. If a real MURR source is later found, pin it in `data/sources.csv` first, with a real
  `url`, `document_version`, `accessed_date`, and `sha256` where applicable — then decide
  separately whether any concrete resulting measurement belongs in one of the
  `data/research/*_candidates.csv` ledgers.
- **Abandon rather than treat a friendly introduction as institutional access.** An intro
  arranged by Caixia Wan to a MURR contact is a conversation opportunity, not a facility
  agreement, an approved protocol, or a biosafety determination.
- **Abandon rather than let this note read as more confident than the ISS precedent this repo
  already audited.** The whole point of tightening the ISS claims in
  `session_claims_2026-08-24_redteam.md` is that an impressive-sounding irradiation
  narrative with no stated dosimeter is exactly the failure mode already caught once. A real
  MURR campaign needs its own dosimeter, its own stated protocol, and its own citation — not
  inherited confidence from either the ISS study or from general knowledge about reactors.
- **If someone proposes writing a `data/research/*_candidates.csv` row from this note alone,
  stop.** There is no dataset, no evidence, and no measurement yet — only a venue lead.

## Next step, if this goes anywhere

After any real conversation with the MURR contact: update this file with what was actually
said (capability, willingness, contact process), and only then evaluate whether it has
become a real, citable `source_id` for `data/sources.csv`.
