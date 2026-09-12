# Note to Caixia Wan — the v1.2 correction

**Status 2026-09-07: sent, twice.** The correction email went to Caixia Wan on **2026-09-05 at 3:33 PM
CDT**, with the text in the next section and the **v1.2 build attached** (tag `preprint-v1.2`, commit
4bf1ba9, 36 pages, 593,010 bytes). A follow-up went on the same thread on **2026-09-06, in the evening
before about 10:50 PM CDT**, with the **v1.2.1 build attached** (tag `preprint-v1.2.1`, commit 3081a58,
37 pages, 593,979 bytes), replacing the file she was told to forward. A separate short message about an
unrelated matter went to her the same evening at about 10:50 PM and carries no preprint content. Both
sends are claim sinks outside this repository: `data/claims_ledger.csv` row `WAN-EMAIL-01` names them,
and any later correction of a number in the sent text must reach her (AGENTS.md, "Correcting a published
number"). The copy she held before was the v1.1 build of 29 August 2026, 08:40 CDT, the working tree of
`d404438`. The 2026-09-06 block that stood here, "still drafted, not sent", was wrong: the mailbox it was
checked against was not the one that sent. No mailbox identifiers or addresses are recorded here.

## Text as sent on 2026-09-05

Hunter's paste from the sent message is the record; the paste ends at the closing line. It is the
copydesk-gated text written against her pinned build (every number read back from the manuscript diff
d404438..3081a58; discipline check all zeros, no hard fails, one accepted advisory) with one
difference from the final draft: the closing line names v1.2, which is what was attached.



```
Hi Ellen,

Thanks for getting back to me, and for asking around about RA positions; I appreciate you looking.

The version you have is v1.1 as built on the morning of 29 August. One bound in section 6.2 was wrong. It counted only one of the two roles a cell can play in the direct radiation term, so its reach is about 1,500 times what your copy states. Two counts, one in Table 4 and one beneath it, that rested on the bound moved with it. The direct term alone decides one accepted move in 206,042, not none, and sixteen moves would reverse without it. The bound was raised in review by Codex on the pull request; the counts follow from re-running the tally with it corrected. The reasoning is shown inline in 6.2, and a test now recomputes the bound from the shipped coefficients and fails on either withdrawn form.

Among the other changes, three matter to a reader. Section 3.10 now says no mechanical term exists in any source file and that a measured modulus could not enter the model as written. Sections 2.5 and 3.8 withdraw the claim that the kNN decision tree operationalizes anything and call the name a coinage. And three cited references were not the papers the text needed and have been replaced, with six uncited entries removed and Malo 2018's authors, title and DOI corrected.

Smaller ones you may notice. Section 2.1 adds Robertson 2012, whose growth result reproduces while its attribution to melanin does not. Section 2.6 cites Turick 2011 as a positive measurement running opposite to the mechanism, quotes Casadevall 2017's own concession, and adds a paragraph on what actually lives in a spent fuel pool. Section 2.2 corrects a continuum-model attribution to Xavier 2005. Section 7.1 and the conclusion now put the melanin term at about ten times the direct radiation term (9.6 in ΔH), where your copy says several orders. And Table 2's phase-locking frequency row is gone; its citations were theory sources and nothing sourced its 0.01 to 1.0 rad per hour range, so the symbol moved to Table 1 with no value assigned. The phenotype boundary in 7.2 and the spent-fuel-pool literature in 7.3 are as you saw them.

I have a short memo on the FeCrAl question if that's useful for the MURR conversation.

v1.2 is attached; it's the one to forward if things move.

Hunter
```

## Follow-up as sent on 2026-09-06

Gated the same way (discipline check all zeros; craft review clean; prose review clean with two style
advisories left unapplied). The date is absolute so this note can quote it in any year.

```
Hi Ellen,

v1.2.1 supersedes the 5 September attachment, with five small corrections and none to the numbers. Please forward this one if things move.

Hunter
```

## What v1.2.1 contains, read from the tagged artifact against the v1.2 build she holds

Corrections, five, none to a reported number: the byline "Hunter Kinder, B.A." (v1.2: "M.A.T.L.");
§6.2's pointer "Software and Data Availability section" (v1.2: "Data and Code Availability"); the Table 4
caption's command printed with two dashes (v1.2: ligated); "the three runs of Table 4" replaced by
"seeds 42, 43 and 44 at 400 MCS (53 603, 68 465 and 83 974 accepted moves respectively; Table 4 shows
seed 42)", which states the per-seed counts for the first time (the phrase survives once in the §8
correction passage); and the abstract's "7.505×10⁻² in ΔH (a 1.5% acceptance bias)". Additions, not
corrections: three Frontiers DOIs on references 51 to 53, and "including as metaphor" dropped from the
Introduction. So "none to the numbers" is true of the build: no number she has was changed; the
per-seed counts and the 1.5% figure were added.

---

## Superseded text (v1.2 form, before the eight-finding review)

Kept as the record; do not send. It failed review on eight points: "v1.1" is nine builds, not
one; "both numbers" named neither; "cannot recur silently" claimed more than the test does;
"automated review" where the manuscript says Codex on the pull request; "four smaller things"
undercounted; "wrong DOIs" understated three replaced papers; the Table 2 row did have citations;
and the attachment line, which in this text meant the untagged build of the day; the sent email also said v1.2, and attached the tagged v1.2 build.

### As it read

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
