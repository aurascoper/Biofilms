#!/usr/bin/env python3
"""Render the consortium panel, the per-MCS frames and the video from seeded frame files.

    python3 tools/render_consortium.py [diagnostics/melanin_ensemble/consortium_seed42] \
        --figures /tmp/biofilms-render/figures --video /tmp/biofilms-render/assets

Reads the frames consortium_frames.jl wrote (occupied sites per arm and MCS, plus a manifest
naming seed, configuration, arms, git commit, Julia version and machine) and writes:

  <figures>/fig_consortium_seed42.{pdf,txt,sha256,png}   a 2 x 3 panel: full and null arms at
                                                          MCS 0, 50 and 100 (the MCS 0 column is
                                                          the shared seeding frame)
  <video>/consortium_seed42_full_vs_null.mp4              every frame to MCS 100, both arms side by side
  <video>/consortium_seed42_full_vs_null.gif              the same, for the README

Voxels are drawn with matplotlib as tools/render_label_trajectory.py draws them, in the same
categorical palette, so the picture needs no ParaView: Julia to make the frames, Python and
ffmpeg to render them. Render from a clean checkout to external output directories.
The image names the committed frame archive and the clean render revision separately.
The .render.json receipt preserves the original run metadata, including its dirty marker;
rendering archived frames does not establish a new simulation origin.
Written 2026-09-24.
"""
import argparse
import csv
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tools"))
from render_label_trajectory import COLOURS, SPECIES  # noqa: E402  the one palette

DEFAULT_FRAMES = REPO / "diagnostics" / "melanin_ensemble" / "consortium_seed42"
NAME = "fig_consortium_seed42"
VIDEO = "consortium_seed42_full_vs_null"
PANEL_MCS = (0, 50, 100)
ARM_TITLE = {"full": "full model", "null": "both radiation-derived terms off"}
PANEL_TITLE = {"full": "full model", "null": "both terms off"}


def manifest(frames):
    m = {}
    for line in (frames / "manifest.txt").read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        m[k] = v
    for key in ("seed", "N", "parcels_per_species", "n_mcs", "arms", "git_commit", "julia_version", "machine"):
        if key not in m:
            raise SystemExit(f"manifest lacks {key}")
    if m.get("seeding_frame_identical_across_arms") != "true":
        raise SystemExit("the manifest does not record a shared seeding frame; refusing to draw the arms side by side")
    return m


def lattice(frames, arm, mcs, n):
    path = frames / f"{arm}_mcs{mcs:03d}.csv.gz"
    if not path.is_file():
        raise SystemExit(f"missing frame {path}")
    a = np.zeros((n, n, n), dtype=np.int16)
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            a[int(r["x"]) - 1, int(r["y"]) - 1, int(r["z"]) - 1] = int(r["species"])
    return a


def draw(ax, species, n, title):
    rgba = np.array([to_rgba(c) for c in COLOURS])
    faces = ax.voxels(species > 0, facecolors=rgba[species], edgecolors=(.13, .16, .2, .38), linewidth=.10, shade=False)
    for face in faces.values():
        face.set_rasterized(True)
    ax.set(xlim=(0, n), ylim=(0, n), zlim=(0, n))
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=25, azim=35)
    ax.set_axis_off()
    ax.set_title(title, fontsize=9)


def clean_render_record(frames):
    """Bind a clean render revision to its archived inputs; never relabel the run."""
    def git(*args):
        return subprocess.check_output(["git", "-C", str(REPO), *args], text=True).strip()

    if git("status", "--porcelain", "--untracked-files=no"):
        raise SystemExit("tracked source files differ from HEAD; commit them before claiming a clean render revision")
    try:
        relative = frames.resolve().relative_to(REPO)
    except ValueError:
        raise SystemExit("frames must be the committed archive inside this checkout")
    tracked = set(git("ls-files", "--", relative.as_posix()).splitlines())
    actual = {p.relative_to(REPO).as_posix() for p in frames.resolve().rglob("*") if p.is_file()}
    if not actual or actual != tracked:
        raise SystemExit("the frame directory must contain exactly its committed files")
    sources = ("tools/render_consortium.py", "tools/render_label_trajectory.py",
               "calibration/tests/test_claims_ledger.py")
    m = manifest(frames)
    return {
        "schema_version": 1,
        "render_git_commit": git("rev-parse", "HEAD"),
        "frame_archive_tree": git("rev-parse", f"HEAD:{relative.as_posix()}"),
        "original_run": {key: m[key] for key in (
            "git_commit", "julia_version", "machine", "seed", "N", "parcels_per_species", "n_mcs", "arms")},
        "source_sha256": {p: hashlib.sha256((REPO / p).read_bytes()).hexdigest() for p in sources},
        "input_sha256": {p: hashlib.sha256((REPO / p).read_bytes()).hexdigest() for p in sorted(actual)},
    }


def external_output(path):
    path = Path(path).resolve()
    if path == REPO or REPO in path.parents:
        raise SystemExit("render outputs must be outside the source checkout; copy them back after validation")
    return path


def provenance(m, record):
    return (f"seed {m['seed']}, N = {m['N']}, {m['parcels_per_species']} parcels per species, to MCS {m['n_mcs']}; "
            f"frame archive {record['frame_archive_tree'][:10]}; render revision {record['render_git_commit'][:10]}")


def panel(frames, m, figdir, record):
    n = int(m["N"])
    arms = m["arms"].split(",")
    if arms != ["full", "null"]:
        raise SystemExit(f"the panel is drawn for arms full,null; the frames declare {arms}")
    plt.rcParams.update({"font.size": 9, "pdf.fonttype": 42, "font.family": "DejaVu Sans"})
    fig = plt.figure(figsize=(10.5, 7.2))
    gs = fig.add_gridspec(2, 3, left=.02, right=.98, top=.86, bottom=.14, hspace=.05, wspace=.02)
    for i, arm in enumerate(arms):
        for j, mcs in enumerate(PANEL_MCS):
            ax = fig.add_subplot(gs[i, j], projection="3d")
            label = "MCS 0, the shared seeding frame" if mcs == 0 else f"MCS {mcs}"
            draw(ax, lattice(frames, arm, mcs, n), n, f"({chr(97 + 3 * i + j)}) {PANEL_TITLE[arm]}, {label}")
    handles = [Patch(facecolor=COLOURS[k + 1], label=name) for k, name in enumerate(SPECIES)]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8, bbox_to_anchor=(0.5, 0.055))
    fig.suptitle("One realisation of the ensemble configuration under two arms of the acceptance rule", fontsize=12, y=0.965)
    fig.text(0.5, 0.045, "Bottom row: both radiation-derived vectors zeroed (beta_ion = 0 and alpha_M = 0); the acceptance rule keeps adhesion and volume only.",
             ha="center", fontsize=8)
    fig.text(0.5, 0.915, provenance(m, record), ha="center", fontsize=8)
    # Correction RAD-FIG-02: the ensemble scope belongs to Figure 1.
    fig.text(0.5, 0.015, "Figure 1 reports the 256-seed ensemble; these frames show seed 42.",
             ha="center", fontsize=8)
    figdir.mkdir(parents=True, exist_ok=True)
    base = figdir / NAME
    fig.savefig(base.with_suffix(".pdf"), dpi=200)
    plt.close(fig)
    subprocess.run(["pdftotext", "-layout", str(base.with_suffix(".pdf")), str(base.with_suffix(".txt"))], check=True)
    txt = base.with_suffix(".txt").read_text(encoding="utf-8").lower()
    if not txt.strip():
        raise SystemExit("pdftotext produced an empty sidecar; refusing to leave it")
    for banned in retracted_figure_terms():
        if banned in txt:
            raise SystemExit(f"the sidecar carries the retracted figure term {banned!r}")
    base.with_suffix(".sha256").write_text(hashlib.sha256(base.with_suffix(".pdf").read_bytes()).hexdigest() + "\n")
    subprocess.run(["pdftoppm", "-png", "-r", "150", "-singlefile", str(base.with_suffix(".pdf")), str(base)], check=True)
    print(f"wrote {base}.pdf/.txt/.sha256/.png")


def retracted_figure_terms():
    import re
    src = (REPO / "calibration" / "tests" / "test_claims_ledger.py").read_text(encoding="utf-8")
    m = re.search(r"RETRACTED_IN_FIGURES\s*=\s*\((.*?)\)", src, re.S)
    if not m:
        raise SystemExit("RETRACTED_IN_FIGURES not found")
    return re.findall(r'"([^"]+)"', m.group(1))


def video(frames, m, videodir, record):
    if shutil.which("ffmpeg") is None:
        raise SystemExit("ffmpeg not found: the video could not be assembled, refusing to leave frames")
    n = int(m["N"])
    arms = m["arms"].split(",")
    n_mcs = int(m["n_mcs"])
    plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans"})
    with tempfile.TemporaryDirectory() as tmp:
        for k in range(n_mcs + 1):
            fig = plt.figure(figsize=(8, 4.4))
            gs = fig.add_gridspec(1, 2, left=.01, right=.99, top=.84, bottom=.04, wspace=.02)
            for i, arm in enumerate(arms):
                ax = fig.add_subplot(gs[0, i], projection="3d")
                draw(ax, lattice(frames, arm, k, n), n, f"{ARM_TITLE[arm]}, MCS {k}")
            fig.suptitle(provenance(m, record), fontsize=8, y=0.96)
            fig.savefig(Path(tmp) / f"frame{k:04d}.png", dpi=110)
            plt.close(fig)
            if k % 20 == 0:
                print(f"frame {k}", file=sys.stderr)
        videodir.mkdir(parents=True, exist_ok=True)
        mp4 = videodir / f"{VIDEO}.mp4"
        gif = videodir / f"{VIDEO}.gif"
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "10", "-i", str(Path(tmp) / "frame%04d.png"),
                        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", str(mp4)], check=True)
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "10", "-i", str(Path(tmp) / "frame%04d.png"),
                        "-vf", "fps=10,scale=720:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse", str(gif)], check=True)
        for p in (mp4, gif):
            if p.stat().st_size == 0:
                raise SystemExit(f"{p} is empty")
            (videodir / (p.name + ".sha256")).write_text(hashlib.sha256(p.read_bytes()).hexdigest() + "\n")
        print(f"wrote {mp4} and {gif}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("frames", nargs="?", default=str(DEFAULT_FRAMES))
    ap.add_argument("--figures", required=True, help="output directory outside the source checkout")
    ap.add_argument("--video", help="video output directory outside the source checkout")
    ap.add_argument("--no-video", action="store_true")
    a = ap.parse_args()
    if not a.no_video and not a.video:
        ap.error("--video is required unless --no-video is given")
    fr = Path(a.frames).resolve()
    figdir = external_output(a.figures)
    videodir = external_output(a.video) if not a.no_video else None
    record = clean_render_record(fr)
    m = manifest(fr)
    panel(fr, m, figdir, record)
    if not a.no_video:
        video(fr, m, videodir, record)
    if clean_render_record(fr) != record:
        raise SystemExit("source or inputs changed during rendering; no render receipt will be written")
    outputs = [(f"preprint/figures/{NAME}.{ext}", figdir / f"{NAME}.{ext}")
               for ext in ("pdf", "txt", "sha256", "png")]
    if not a.no_video:
        outputs.extend((f"assets/{VIDEO}.{ext}", videodir / f"{VIDEO}.{ext}")
                       for ext in ("mp4", "gif", "mp4.sha256", "gif.sha256"))
    record["output_sha256"] = {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in outputs}
    receipt = figdir / f"{NAME}.render.json"
    receipt.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {receipt}")
