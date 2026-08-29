# Negative-control fixtures

Each file here is a **known-bad input** that a guard in
`calibration/tests/test_claims_ledger.py` must detect. A guard that finds nothing
looks identical to a guard that can find nothing (AGENTS.md rule 1); these are what
make the difference visible.

They are committed rather than recovered with `git show <sha>:<path>`. That idiom
was used three times — `5980dc5`, `9319d43`, `e24dbec` — and each time the same
exposure was written down and left in place: a squash-merge makes the pinned commit
unreachable, and the control silently degrades into a skip, which rule 2 says to read
as uncovered surface. A committed file cannot become unreachable.

| fixture | taken from | the guard it must fail |
|---|---|---|
| `modeling_radiotrophic_fitness_prerevision.md` | `5980dc5:preprint/modeling_radiotrophic_fitness.md` | >= 18 `delete` claims must be detected in the pre-revision manuscript |
| `wan_meeting_handout_prefix.tex` | `9319d43:preprint/wan_meeting_handout.tex` | `HANDOUT-02` must be detected |
| `fig1_radial_stratification_prefix.pdf` | `e24dbec:preprint/figures/fig1_radial_stratification.pdf` | `radiotrophic niche` / `radiosensitive core` must be detected |
| `fig2_melanin_accumulation_prefix.pdf` | `e24dbec:preprint/figures/fig2_melanin_accumulation.pdf` | `... are radiotrophic (melanin-mediated energy gain)` must be detected |

Do not regenerate or "fix" these files. Their value is that they are wrong.
