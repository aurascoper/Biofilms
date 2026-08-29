# Draft note to Caixia Wan — v1.1 figure correction

Not sent. Send it before Andrew Hoffman opens the version she has.

---

**Subject:** Corrected version of the preprint (v1.1) — figure error in the copy you have

Dr. Wan,

Thank you again for Thursday, and for offering to send the preprint to Andrew Hoffman.

Before you do: I found an error in the figures of the version you have. Attached is v1.1,
which fixes it.

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

No number changed. The underlying simulation is the same run, and the values quoted in the
captions reproduce exactly. What changed is five pieces of text inside two images, plus a
correction note in the Software and Data Availability section recording it.

Two other additions in v1.1, both prompted by the gap you identified:

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

I have also written a short memo for Dr. Hoffman, if it is useful to send alongside: one page
on what the model computes, one on why a materials scientist rather than a microbiologist,
and one page that is entirely questions, since he knows what MURR can do and I do not.

Thank you again — for the time, and for the assistantship conversation.

Hunter Kinder
