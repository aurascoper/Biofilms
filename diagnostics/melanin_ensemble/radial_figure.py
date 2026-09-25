#!/usr/bin/env python3
"""Render the radial ensemble as the manuscript's Figure 1, from the committed CSV alone.

    python3 diagnostics/melanin_ensemble/radial_figure.py [radial.csv] [--out preprint/figures]

WHAT IT SHOWS. Panel (a): per-species mean r/R at MCS 100 under the four arms, as the
ensemble mean with a 95% interval of that mean over 256 seeds, beside the seeding frame at
MCS 0. Panel (b): the paired per-seed difference full - no_mel for C. neoformans and
C. sphaerospermum, the two species the pre-registered test P1 named, with the minimum effect
of interest drawn; and full - no_rad for every species with the equivalence bound drawn.
Panel (c): the C. neoformans - B. subtilis separation per seed under the full arm, sorted,
with seed 42 marked, because that one draw was Figure 1 of versions 1.0 to 1.2.

WHY MATPLOTLIB. The other figure renderers in tools/ use it, and it needs no Julia
package server: the CSV is the only input. Written 2026-09-24.

THE SIDECARS ARE WRITTEN HERE, in the order the guards read them: the PDF, its pdftotext
-layout text, its sha256 (of the PDF), and the PNG. A refusal on any of them leaves no
figure behind a stale sidecar.
"""
import argparse
import csv
import hashlib
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[2]
DEFAULT_CSV = REPO / "diagnostics" / "melanin_ensemble" / "radial_n40_p6_seeds42-297.csv"
NAME = "fig1_radial_ensemble"
ARMS = ["full", "no_rad", "no_mel", "null"]
ARM_LABEL = {"full": "full", "no_rad": r"$\beta_{ion}=0$", "no_mel": r"$\alpha_M=0$", "null": "both off"}
SPECIES = [(1, "C. neoformans"), (2, "D. radiodurans"), (3, "C. sphaerospermum"), (4, "B. subtilis"),
           (5, "A. niger"), (6, "S. oneidensis"), (7, "O. intermedium AM7")]
MEOI = 0.05          # the pre-registered minimum effect and equivalence bound, r/R
PUBLISHED_SEED = 42
# Validated with the dataviz palette validator on a light surface (adjacent CVD dE >= 8).
COLOURS = {"full": "#2a78d6", "no_rad": "#eb6834", "no_mel": "#3aa655", "null": "#6b6b68"}
INK = "#2f2f2e"


def retracted_figure_terms():
    """The figure vocabulary the guards refuse, read from its one source so this file does not
    carry the words itself; tests/manuscript_claims_tests.jl reads it the same way."""
    src = (REPO / "calibration" / "tests" / "test_claims_ledger.py").read_text(encoding="utf-8")
    m = re.search(r"RETRACTED_IN_FIGURES\s*=\s*\((.*?)\)", src, re.S)
    if not m:
        raise SystemExit("RETRACTED_IN_FIGURES not found in test_claims_ledger.py")
    return [t for t in re.findall(r'"([^"]+)"', m.group(1))]


def read(path):
    """species -> seed -> arm -> r/R at the read MCS, plus the seeding frame and header."""
    meta, at, frame = {}, defaultdict(lambda: defaultdict(dict)), defaultdict(dict)
    with open(path, encoding="utf-8") as fh:
        rows = []
        for line in fh:
            if line.startswith("#"):
                for k, v in (t.split("=", 1) for t in line[1:].split() if "=" in t):
                    meta.setdefault(k, v)
                continue
            rows.append(line)
    read_at = int(meta["read_at_mcs"])
    for r in csv.DictReader(rows):
        sp, seed, arm, mcs = int(r["species"]), int(r["seed"]), r["arm"], int(r["mcs"])
        v = float(r["mean_r_over_R"])
        if mcs == read_at:
            if arm in at[sp][seed]:
                raise SystemExit(f"duplicate row: seed {seed} arm {arm} species {sp}")
            at[sp][seed][arm] = v
        elif mcs == 0:
            if seed in frame[sp] and frame[sp][seed] != v:
                raise SystemExit(f"seed {seed} species {sp}: the seeding frame differs across arms")
            frame[sp][seed] = v
    arms = meta["arms"].split(",")
    if arms != ARMS:
        raise SystemExit(f"{path} declares arms {arms}; this figure is drawn for {ARMS}")
    for sp, _ in SPECIES:
        for seed, per_arm in at[sp].items():
            missing = [a for a in ARMS if a not in per_arm]
            if missing:
                raise SystemExit(f"seed {seed} species {sp} lacks {missing}")
    return meta, at, frame, read_at


def ci95(x):
    x = np.asarray(x)
    return 1.96 * x.std(ddof=1) / np.sqrt(len(x))


def render(csv_path, outdir):
    try:
        subprocess.run(["pdftotext", "-v"], capture_output=True, check=False)
    except FileNotFoundError:
        raise SystemExit("pdftotext not found: the .txt sidecar could not be written, refusing to render")
    meta, at, frame, read_at = read(csv_path)
    seeds = sorted(at[1])
    n = len(seeds)
    if PUBLISHED_SEED not in seeds:
        raise SystemExit(f"seed {PUBLISHED_SEED} is not in the file; the marker would lie")

    plt.rcParams.update({"font.size": 9, "pdf.fonttype": 42, "font.family": "DejaVu Sans"})
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 5.0), gridspec_kw={"width_ratios": [1.45, 1.15, 1.0]})
    fig.subplots_adjust(left=0.12, right=0.985, top=0.80, bottom=0.34, wspace=0.55)

    # (a) per-species mean r/R per arm, with the seeding frame
    ax = axes[0]
    y = np.arange(len(SPECIES))
    offsets = np.linspace(-0.3, 0.3, len(ARMS) + 1)
    for j, (sp, name) in enumerate(SPECIES):
        f = np.array([frame[sp][s] for s in seeds])
        ax.errorbar(f.mean(), y[j] + offsets[0], xerr=ci95(f), fmt="s", color=INK, ms=3.5, capsize=2,
                    label="seeding frame, MCS 0" if j == 0 else None)
        for k, arm in enumerate(ARMS):
            v = np.array([at[sp][s][arm] for s in seeds])
            ax.errorbar(v.mean(), y[j] + offsets[k + 1], xerr=ci95(v), fmt="o", color=COLOURS[arm], ms=3.5,
                        capsize=2, label=f"{ARM_LABEL[arm]}, MCS {read_at}" if j == 0 else None)
    ax.set_yticks(y)
    ax.set_yticklabels([n_ for _, n_ in SPECIES], fontstyle="italic")
    ax.invert_yaxis()
    ax.set_xlabel(f"mean r/R over parcels\n(ensemble mean, 95% interval of the mean, {n} seeds)")
    ax.set_title(f"(a) Every arm relaxes inward together")
    ax.axvline(2 / 3 * 0.8, color=INK, lw=0.6, ls=":", label="centre-only null, 2/3 of 0.8R")
    ax.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.30), fontsize=7.5, ncol=2, frameon=False)
    ax.grid(axis="x", lw=0.3, alpha=0.5)

    # (b) the two primary tests
    ax = axes[1]
    rows = []
    for sp, name in SPECIES:
        d = np.array([at[sp][s]["full"] - at[sp][s]["no_rad"] for s in seeds])
        rows.append((f"{name}\nfull $-$ ($\\beta_{{ion}}=0$)", d, COLOURS["no_rad"]))
    for sp, name in [(1, "C. neoformans"), (3, "C. sphaerospermum")]:
        d = np.array([at[sp][s]["full"] - at[sp][s]["no_mel"] for s in seeds])
        rows.append((f"{name}\nfull $-$ ($\\alpha_M=0$), P1", d, COLOURS["no_mel"]))
    yy = np.arange(len(rows))
    for j, (lab, d, c) in enumerate(rows):
        half = 1.6448536 * d.std(ddof=1) / np.sqrt(len(d))
        ax.errorbar(d.mean(), yy[j], xerr=half, fmt="o", color=c, ms=3.5, capsize=2)
    ax.axvspan(-MEOI, MEOI, color=INK, alpha=0.06, lw=0)
    ax.axvline(0, color=INK, lw=0.6)
    ax.set_yticks(yy)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7.2, fontstyle="italic")
    ax.invert_yaxis()
    ax.set_xlim(-MEOI * 1.4, MEOI * 1.4)
    ax.set_xlabel("paired mean difference in r/R, 90% interval")
    ax.set_title("(b) Switching a term off moves nothing")
    ax.text(0.5, 0.995, "shaded: the pre-registered bound, ±0.05 r/R", transform=ax.transAxes, ha="center", va="top", fontsize=7, color=INK)
    ax.grid(axis="x", lw=0.3, alpha=0.5)

    # (c) CN - BS separation per seed, full arm, sorted, seed 42 marked
    ax = axes[2]
    sep = {s: at[1][s]["full"] - at[4][s]["full"] for s in seeds}
    order = sorted(seeds, key=lambda s: sep[s])
    vals = np.array([sep[s] for s in order])
    ax.bar(np.arange(n), vals, color=COLOURS["full"], width=1.0, lw=0)
    k = order.index(PUBLISHED_SEED)
    ax.plot([k], [vals[k]], marker="v", color="#eb6834", ms=7)
    ax.annotate(f"seed {PUBLISHED_SEED}: {sep[PUBLISHED_SEED]:+.2f}", xy=(k, vals[k]), xytext=(-70, 12),
                textcoords="offset points", fontsize=7.5, color="#eb6834",
                arrowprops=dict(arrowstyle="-", color="#eb6834", lw=0.6))
    ax.axhline(0, color=INK, lw=0.6)
    ax.set_xlabel(f"seeds, sorted ({n})")
    ax.set_ylabel("C. neoformans $-$ B. subtilis (r/R), full arm")
    ax.set_title("(c) Seed 42 was one draw")
    ax.text(0.02, 0.96, f"mean {vals.mean():+.3f}, sd {vals.std(ddof=1):.2f}", transform=ax.transAxes,
            va="top", fontsize=7.5)
    ax.set_xticks([])

    fig.suptitle("Neither radiation-derived term moves a parcel: four arms on the same 256 seeds", fontsize=12, y=0.975)
    fig.text(0.5, 0.905, f"N = 40, six parcels per species, read at MCS {read_at}, seeds {seeds[0]} to {seeds[-1]}; "
             f"paired within seed; produced by julia {meta.get('julia_version','?')} on {meta.get('machine','?')}",
             ha="center", fontsize=8)
    fig.text(0.5, 0.035, "Lattice units divided by R = N/2. Arms zero CPMParams vectors only. A model output about the acceptance rule, "
             "with no organism in it.", ha="center", fontsize=8)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    base = outdir / NAME
    fig.savefig(base.with_suffix(".pdf"), dpi=220)
    plt.close(fig)
    subprocess.run(["pdftotext", "-layout", str(base.with_suffix(".pdf")), str(base.with_suffix(".txt"))], check=True)
    txt = base.with_suffix(".txt").read_text(encoding="utf-8")
    if not txt.strip():
        raise SystemExit("pdftotext produced an empty sidecar; refusing to leave it")
    for banned in retracted_figure_terms():
        if banned in txt.lower():
            raise SystemExit(f"the sidecar carries the retracted figure term {banned!r}")
    digest = hashlib.sha256(base.with_suffix(".pdf").read_bytes()).hexdigest()
    base.with_suffix(".sha256").write_text(digest + "\n")
    # PNG last, rasterised from the PDF so the two cannot differ.
    subprocess.run(["pdftoppm", "-png", "-r", "160", "-singlefile", str(base.with_suffix(".pdf")), str(base)], check=True)
    print(f"wrote {base}.pdf/.txt/.sha256/.png  sha256 {digest}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="?", default=str(DEFAULT_CSV))
    ap.add_argument("--out", default=str(REPO / "preprint" / "figures"))
    a = ap.parse_args()
    render(a.csv, a.out)
