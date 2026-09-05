# Note to Caixia Wan — the v1.2 correction

**Status 2026-09-05: drafted in Gmail, not sent.** To Caixia Wan, subject "Preprint v1.2,
correcting the 29 August v1.1 you have". The copy she holds is the v1.1 build of 29 August 2026,
08:40 CDT, the working tree of `d404438` (bound −5×10⁻⁵, Table 4 with "none of 206 042",
"several orders" in §7.1, figures already corrected). The attachment is the CI build of the
commit tagged `preprint-v1.2`, added by hand because no tool here attaches a file of that size.
When it goes, replace this block with the sent date: a sent email is a claim sink and belongs in
the dependents enumeration (AGENTS.md, "Correcting a published number"). No mailbox identifiers
or addresses are recorded here; the mailbox holds them.

## Final text (v1.2 form)

Everything in it was checked against the diff f72aabb..26b3a14 of the manuscript: the
withdrawn bound lived in section 6.2 of v1.1, not the abstract (the v1.1 abstract carried no
figure); sections 7.2 and 7.3 have no content change; the four listed changes are the ones
outside 6.2 a reader would notice. Gated through copydesk (no hard fails) and prose-craft-2.

```
Hi Ellen,

Thanks for getting back to me, and for asking around about RA positions; I appreciate you looking.

The version you have is v1.1, and I've since corrected two numbers in it. The direct radiation term's effect on acceptance was understated. The bound I'd published counted only one of the two roles a cell can play in the energy term, so the real reach is about 1,500 times larger than the bound stated in section 6.2 of your copy. Both numbers were caught in automated review of the pull request and are corrected in section 6.2 with the reasoning shown, and there is now a test that recomputes the bound from the code, so that class of error cannot recur silently.

Four smaller things moved with it, and I'd rather you heard them from me than found them. Section 2.1 adds Robertson 2012, whose growth result reproduces while its attribution to melanin does not. Section 2.6 adds a paragraph on what actually lives in a spent fuel pool. Table 2 lost its phase-locking frequency row, which had no literature source. And several reference entries with wrong DOIs were corrected. The phenotype boundary in 7.2 and the spent-fuel-pool literature in 7.3 are as you saw them.

I have a short memo on the FeCrAl question if that's useful for the MURR conversation.

v1.2 is attached; it's the one to forward if things move.

Hunter
```

---

## Earlier draft (v1.1 form, superseded by the text above)

Superseded twice since drafting, and the second time it said something false. It
was written for v1.1, which corrected two figures and no numbers. v1.2 corrects numbers,
including the radiation count this draft quotes. Kept as the record of what the note said
before the diff was checked; do not send any part of it.

---

**Subject:** Corrected version of the preprint (v1.2) — figure error, and a number, in the copy you have

Dr. Wan,

Thank you again for Thursday, and for offering to send the preprint to Andrew Hoffman.

Before you do: I found an error in the figures of the version you have, and then, reviewing the
section that error sat next to, an error in one of its numbers. Attached is v1.2, which fixes
both. The figure error is the one I would have written to you about on its own; the numerical one
is more serious and I would rather you heard it from me than found it.

Figure 2 in v1.0 carries a line inside the plot image reading "C. neoformans,
C. sphaerospermum are radiotrophic (melanin-mediated energy gain)". That contradicts the
paper's own Section 2.6, which states that radiotrophy is not established for any of the
seven species modelled, and it contradicts the caption directly beneath it. Figure 1 has a
shaded band labelled "radiotrophic niche", on the wrong side of the plot as well.

What happened is that I corrected the figure-generating code two weeks ago and never
regenerated the committed images, and nothing in the test suite could open a figure to
notice the difference. The prose was audited to the ground; the text baked into the images
was not. It is the paper's own argument happening to the paper, which is an uncomfortable
thing to have to write and the reason I would rather you had the corrected copy.

No number changed *in the figures*. The underlying simulation is the same run, and the values
quoted in the captions reproduce exactly. What changed is five pieces of text inside the two
images — both plot titles, the two shaded-band labels in Figure 1, and the annotation quoted
above — plus a correction note in the Software and Data Availability section recording it.

Numbers did change in Section 6.2, and that is the second correction. An automated review of the
code found that the bound I used to argue one of that section's results was wrong by a factor of
about 1500 — I had enumerated only one of the two roles a cell can play in the energy term, and
then, on a second pass, had bounded the corrected quantity over species when it is an extremum
over pairs of them. Two published counts moved with it. The section now states what the
measurement supports rather than what the bad bound predicted, and the code carries a test that
recomputes the bound and fails if the paper and the coefficients disagree. The corrections are
recorded in the Data and Code Availability section and in the claims ledger.

A few additions have gone in alongside the correction, most of them prompted by the gap you
identified:

- **Section 7.3** now discusses biofilms in nuclear facilities — spent fuel pools and cooling
  circuits — alongside the environmental remediation framing. The argument is about regime:
  the paper already notes that reactor irradiation and contaminated-site dose rates are about
  ten orders of magnitude apart, and a facility sits between them, which is where the model's
  inputs are actually obtainable. There is published work on biofilms retaining Co-60 on
  stainless steel and titanium coupons in spent fuel pools, and one of the isolates in that
  literature is *B. subtilis*, already one of the seven species in the model.

- **The Ethics statement** now carries the taxonomy of *Ochrobactrum intermedium* AM7, which
  was reclassified into *Brucella* in 2020 and is served by NCBI as *Brucella intermedia*.
  Both names are on the record now. No containment determination is claimed — that follows
  the strain and belongs to an institution — but anyone reading the paper should meet the
  name change in the paper rather than in a biosafety committee.

- **Section 6.2** has a new table. It asks, for every accepted move in the simulation, whether
  removing one term of the model's energy function would have reversed it. Removing the direct
  radiation term reverses sixteen of 206,042 accepted moves across three seeds, and in fifteen of
  those an adhesion or volume term was independently decisive as well, so the term is the sole
  decider of exactly one. Two-thirds to three-quarters of moves are reversed by removing no single
  term at all, so the dynamics are carried by the sum rather than by any one component.

  *(This bullet said "reversed none of 206,042" when the note was drafted for v1.1. That was the
  withdrawn count, and it is the number the v1.2 correction moves. The distinction that matters is
  absorption rather than absence: the term reaches far enough to have decided sixteen moves, and
  is usually not the only term that could have.)*

There is also some tightening of terminology — one coined term is no longer used before the
paragraph that explains it was coined — which changes no result.

**Two separate things, so you do not have to work out which is which.** One is a measurement in
your own lab: the unirradiated biosorption assay in the handout, ICP readout, no reactor time and no
scheduling. The other is forwarding the memo below to Dr. Hoffman, which asks him about a material
and a facility and needs nothing from you but the introduction. They are independent — either can go
ahead without the other, and neither is a precondition for the other.

I have also written a memo for Dr. Hoffman, if it is useful to send alongside. Four short
landscape pages: what the model computes and what it does not, the published work on biofilms
in spent fuel pools, why I am writing to a materials scientist at all — his published corrosion
work qualifies alloys against steam, hydrothermal chemistry and hydrogen permeation, and a
biofilm is none of those — and a final page that is entirely questions, since he knows what
MURR can do and I do not.

Thank you again — for the time, and for the assistantship conversation.

Hunter Kinder
