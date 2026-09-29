"""Check the 2026-09-29 prose corrections against the stored radial ensemble.

Scope: the named sentences, not an inference that all manuscript claims are
correct. The control contains verbatim excerpts extracted from 24cd3ec; it
must fail for all seven defects before the corrected manuscript can pass.
No simulation is rerun. The interval follows radial_analyse.jl's normal 90% CI.
"""
from __future__ import annotations

import csv
import json
import math
import re
import statistics
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TEX = REPO / "preprint/modeling_radioresistance_and_radiotropic_fitness.tex"
CSV = REPO / "diagnostics/melanin_ensemble/radial_n40_p6_seeds42-297.csv"
CONTROL = Path(__file__).parent / "fixtures/presend_24cd3ec_excerpts.json"


def _between(source, start, end):
    assert source.count(start) == 1, f"missing or ambiguous landmark: {start}"
    tail = source.split(start, 1)[1]
    assert end in tail, f"missing closing landmark: {end}"
    return start + tail.split(end, 1)[0]


def manuscript_blocks(source):
    paragraphs = source.split("\n\n")
    dose = [p for p in paragraphs if r"\mathrm{min}^{-1}" in p
            and r"\mathrm{yr}^{-1}" in p]
    return {
        "abstract": _between(source, r"\begin{abstract}", r"\end{abstract}"),
        "result": _between(source, r"\paragraph{Result.}",
                           r"\paragraph{What it establishes.}"),
        "dose": "\n\n".join(dose),
        "authors": _between(source, r"\section*{Author Contributions}",
                            r"\section*{Software and Data Availability}"),
        "caption": _between(source, r"\includegraphics[width=0.99\linewidth]{fig_consortium_seed42}",
                            r"\label{fig:consortium}"),
    }


def ensemble_facts():
    with CSV.open(encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(l for l in fh if not l.startswith("#"))
                if r["mcs"] == "100"]
    radii = {(int(r["seed"]), r["arm"], int(r["species"])):
             float(r["mean_r_over_R"]) for r in rows}
    seeds, arms = range(42, 298), ("full", "no_rad", "no_mel", "null")
    expected = {(s, a, sp) for s in seeds for a in arms for sp in range(1, 8)}
    assert len(rows) == len(radii) and set(radii) == expected, "incomplete or duplicate pairing"
    means, intervals = {}, {}
    for sp in range(1, 8):
        differences = [radii[s, "full", sp] - radii[s, "no_rad", sp] for s in seeds]
        means[sp] = statistics.mean(differences)
        half = statistics.NormalDist().inv_cdf(.95) * statistics.stdev(differences) / math.sqrt(len(differences))
        intervals[sp] = (means[sp] - half, means[sp] + half)
    separation = [statistics.mean(radii[s, a, 1] - radii[s, a, 4] for s in seeds) for a in arms]
    melanin = [statistics.mean(radii[s, "full", sp] - radii[s, "no_mel", sp]
                              for s in seeds) for sp in (1, 3)]
    return means, intervals, separation, melanin


def correction_errors(blocks, facts):
    """Named defects in the stated corrections; missing expected text also fails."""
    means, intervals, separation, melanin = facts
    text = {k: " ".join(v.split()) for k, v in blocks.items()}
    result, errors = text["result"], set()
    farthest = max(intervals, key=lambda sp: max(map(abs, intervals[sp])))
    interval = re.search(r"farthest from zero.*?\$\[(-?[0-9.]+),\s*(-?[0-9.]+)\]\$"
                         r".*?S\.~oneidensis", result)
    if not (farthest == 6 and interval and
            tuple(map(float, interval.groups())) == tuple(round(x, 4) for x in intervals[6])):
        errors.add("interval_rank")
    ratio = .05 / max(map(abs, means.values()))
    bound = re.search(r"the bound is about thirty times ([^.]+)", result)
    if not (round(ratio / 10) * 10 == 30 and bound and "largest" in bound[1]
            and "mean" in bound[1]):
        errors.add("bound_ratio")
    sep = re.search(r"separation\s+(.*?)(?:standard deviation|over seeds)", result)
    endpoints = re.findall(r"-0\.\d+", sep[1]) if sep else []
    if list(map(float, endpoints[:2])) != [round(min(separation), 4), round(max(separation), 4)]:
        errors.add("separation_range")
    abstract = text["abstract"]
    contrast = re.search(r"full(?:\s+|-)+minus(.{0,180})", abstract, re.I)
    if not (contrast and ("melanin" in contrast[1] or r"\alpha_M" in contrast[1])
            and all(f"{x:+.4f}" in abstract for x in melanin)):
        errors.add("contrast_direction")
    # 1 Gy/min divided by 1 mGy/year, with a 365.25-day year.
    orders = math.log10(1000 * 60 * 24 * 365.25)
    dose_paragraphs = blocks["dose"].split("\n\n")
    if not (8.5 < orders < 9 and len(dose_paragraphs) == 2 and
            all("nearly nine orders of magnitude" in " ".join(p.split())
                for p in dose_paragraphs)):
        errors.add("dose_rate_orders")
    if re.search(r"both authors approved|authors have approved", text["authors"], re.I):
        errors.add("author_approval")
    caption = text["caption"]
    if ("Every reported number in this paper" in caption or
            not re.search(r"Figure~\\ref\{fig:radial\} reports the "
                          r"256-seed (?:radial )?ensemble", caption)):
        errors.add("caption_scope")
    return errors


def test_presend_corrections_match_the_stored_ensemble_and_unit_conversion():
    assert not correction_errors(manuscript_blocks(TEX.read_text(encoding="utf-8")), ensemble_facts())


def test_the_actual_uncorrected_excerpts_fail_for_all_seven_reasons():
    control = json.loads(CONTROL.read_text(encoding="utf-8"))
    assert control["revision"] == "24cd3ec"
    assert correction_errors(control["excerpts"], ensemble_facts()) == {
        "interval_rank", "bound_ratio", "separation_range", "contrast_direction",
        "dose_rate_orders", "author_approval", "caption_scope",
    }
