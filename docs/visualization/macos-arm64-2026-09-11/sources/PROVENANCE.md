# Evaluated decay data for Lu-177 — source audit, closed

`Lu-177.lara.txt` is the LNHB/DDEP evaluated nuclide table, retrieved 2026-09-11 from
`http://www.lnhb.fr/nuclides/Lu-177.lara.txt` (the site does not serve HTTPS).

```
sha256  5e68ccb03fedb7087677594bcc2ecc15d89fb285947f314bad4c467af27e45a8
```

**The file has CRLF terminators and `sources/.gitattributes` sets `* -text` to preserve
them.** Without it git normalised CRLF to LF on commit, so the stored blob hashed to
`09c78fdf...` while this document and `decay_reference.py` both said `5e68ccb0...` — the
working copy matched the pin and **every clone did not**. A content pin the VCS silently
rewrites is worse than no pin, because it fails only for the people who check it. Verify with
`git show HEAD:<path> | shasum -a 256`, not just against the working copy.

The companion `Lu-177_com.pdf` (Comments on evaluation, sha256
`4b7d9bc8364bad00b9f3f1d12d85cab1dfb979a99cb93952b54d2a7e525c6c0f`) is **not** committed —
it is 385 KB and adds attribution rather than data.

## The pinned values, verbatim from the table

```
Nuclide ; Lu-177
Daughter(s) ; (B-) ; Hf-177 ; 100
Q- ; 496.8
Half-life (d) ; 6.6443 ; 0.0009
Half-life (s) ; 574.07E3 ; 0.08E3
Decay constant (1/s) ; 1.20743E-6 ; 0.00016E-6
Specific activity (Bq/g) ; 4.1081E15 ; 0.0006E15
Reference ; CEA/LNE-LNHB - 2025
```

**Version: `CEA/LNE-LNHB - 2025`.** Evaluators M.A. Kellett and X. Mougeot, Université
Paris-Saclay, CEA, List, Laboratoire National Henri Becquerel (LNE-LNHB), F-91120 Palaiseau.
The comments record an initial evaluation in February 2002 by F.G. Kondev (ANL) and a complete
re-evaluation in April 2015 by M.A. Kellett (2016Ke06).

## What this settles

`decay_reference.py` used **6.6443 d** and **λ = 1.20743e-6 /s**, chosen because they were the
self-consistent pair — the accepted plan carried 6.6443 d in one section and **6.647 d** in its
decay equation, and only the first reproduces the stated λ. That reasoning was arithmetic, and
the code said so, carrying a standing confirmation requirement.

**Both are the evaluated values.** The table states them directly, with uncertainties the code
did not have: ±0.0009 d and ±0.00016e-6 /s, a relative uncertainty of 1.35e-04. **6.647 d is
not the evaluated value.**

Internal consistency of the table itself:

| | |
|---|---|
| `6.6443 d x 86400` | 574067.52 s against a stated 574070 s |
| `ln2 / 574070 s` | 1.207426e-06 /s against a stated 1.20743e-06 |
| `ln2 / 1.20743e-6 / 86400` | 6.64431 d against a stated 6.6443 d |

The 574067.52 vs 574070 gap is the table's own rounding of the seconds figure to six
significant figures; it is well inside the stated ±80 s.

## Two further confirmations, not sought but obtained

The plan's Phase 1b constants, previously "pinned from memory" and never verified:

| plan | table |
|---|---|
| γ 112.95 keV (6.22%) | 112.95005 keV, 6.223 ± 0.032 % |
| γ 208.37 keV (10.43%) | 208.3661 keV, 10.425 ± 0.035 % |
| β⁻ Emax 497 keV | `Q- ; 496.8` |
| daughter stable | `Daughter(s) ; (B-) ; Hf-177 ; 100` — Hf-177 is stable |

## What this does NOT settle

A sourced decay constant does not give this repository an isotope identity. The source-term
gate stands, the material path remains elemental rather than isotopic, and
`decay_reference.py` still computes **no binding and no dose**. Naming the isotope is
legitimate for this constant and nothing else; the files, fields and receipt continue to say
*decay reference*.

`a0(x)` remains a declared hypothetical input. Nothing here makes it a measurement.
