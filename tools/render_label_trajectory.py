#!/usr/bin/env python3
"""Render Figure 5 and stage the manuscript from a verified labelled run.

Uses the binary indicators/identity counts produced by Julia, never c or s.
The 3D panels are CPU raster projections of actual voxel faces with one camera.
No ray tracer or smoothing/interpolation of categorical labels is involved.
Requires numpy, matplotlib, h5py and the repository's existing LaTeX/Poppler tools.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = "modeling_radioresistance_and_radiotropic_fitness.tex"
FIGURE = "fig5_label_trajectory"
COLOURS = ["#ffffff", "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6"]
SPECIES = ["C. neoformans", "D. radiodurans", "C. sphaerospermum", "B. subtilis",
           "A. niger", "S. oneidensis", "O. intermedium"]


def sha(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def verified_manifest(run: Path) -> dict:
    manifest_path = run / "run_manifest.json"
    if sha(manifest_path) != (run / "run_manifest.sha256").read_text().strip():
        raise ValueError("manifest receipt mismatch")
    manifest = json.loads(manifest_path.read_text())
    if manifest["status"] not in {"trajectory_verified", "manuscript_built"}:
        raise ValueError("complete trajectory postflight is required before rendering")
    if not manifest["verification"]["determinism_executed"]:
        raise ValueError("determinism verification was not executed")
    for rel, expected in manifest["artifacts"].items():
        if sha(run / rel) != expected:
            raise ValueError(f"artifact hash mismatch: {rel}")
    for rel, expected in manifest["figure1_hashes"].items():
        if sha(ROOT / rel) != expected:
            raise ValueError("Figure 1 changed")
    return manifest


def save_manifest(run: Path, manifest: dict) -> None:
    manifest["artifacts"] = {str(p.relative_to(run)): sha(p)
                             for p in sorted(run.rglob("*")) if p.is_file()
                             and p.relative_to(run).as_posix() not in {"run_manifest.json", "run_manifest.sha256"}}
    (run / "run_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    (run / "run_manifest.sha256").write_text(sha(run / "run_manifest.json") + "\n")


def xyz(dataset) -> np.ndarray:
    # HDF5.jl declares that h5py sees zyx. Conversion is explicit and is also
    # checked against asymmetric producer probes for the snapshot panels.
    return np.asarray(dataset).transpose(2, 1, 0)


def snapshot(run: Path, mcs: int):
    with h5py.File(run / "snapshots" / f"snap_mcs{mcs:06d}.h5") as f:
        order = f.attrs["dataset_axis_order_h5py"]
        if isinstance(order, bytes):
            order = order.decode()
        if order != "zyx":
            raise ValueError("unknown HDF5 axis order")
        if int(f.attrs["mcs"]) != mcs:
            raise ValueError("snapshot MCS mismatch")
        ids = xyz(f["lattice/cell_id"])
        for x, y, z, expected in np.asarray(f["orientation_probes"]).T:
            if ids[x, y, z] != expected:
                raise ValueError("orientation probe mismatch")
        return xyz(f["lattice/species_id"]), int(f.attrs["live_registry_count"])


def render(run: Path, manifest: dict) -> Path:
    output = run / "figures"
    output.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13,
                         "axes.titlesize": 15, "axes.labelsize": 12,
                         "pdf.fonttype": 42})
    fig = plt.figure(figsize=(10.8, 8.1), facecolor="white")
    grid = fig.add_gridspec(2, 6, left=.06, right=.94, top=.89, bottom=.23,
                           hspace=.37, wspace=.8, height_ratios=(1, 1.1))
    rgba = np.array([to_rgba(c) for c in COLOURS])
    for k, mcs in enumerate((0, 30, 100)):
        species, count = snapshot(run, mcs)
        ax = fig.add_subplot(grid[0, 2*k:2*k+2], projection="3d")
        faces = ax.voxels(species > 0, facecolors=rgba[species],
                          edgecolors=(.13, .16, .2, .38), linewidth=.10, shade=False)
        for artist in faces.values():
            artist.set_rasterized(True)
        ax.set(xlim=(0, 40), ylim=(0, 40), zlim=(0, 40))
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=25, azim=35)
        ax.set_axis_off()
        ax.set_title(f"({chr(97+k)})  MCS {mcs}", pad=0)
        ax.text2D(.5, -.02, f"{count} live parcels", ha="center", transform=ax.transAxes, fontsize=12)
    with h5py.File(run / "analysis" / "label_dynamics.h5") as f:
        mask = xyz(f["interior_mask"]).astype(bool)
        # N=40 has two central slices; declare the lower one, z index 20 in Julia.
        z = 19
        for k, (field, title) in enumerate((("transitions", "(d)  Parcel-ID transitions"),
                                           ("returns", "(e)  Returns to initial parcel ID"))):
            values = xyz(f[f"cadence_1/parcel/{field}"])
            plane = np.ma.array(values[:, :, z].T, mask=~mask[:, :, z].T)
            ax = fig.add_subplot(grid[1, k*3:k*3+3])
            cmap = plt.colormaps["YlGnBu"].copy()
            cmap.set_bad("#e4e9ed")
            image = ax.imshow(plane, origin="lower", extent=(0, 40, 0, 40), interpolation="nearest",
                              cmap=cmap, vmin=0, vmax=max(1, int(plane.max())))
            ax.set_title(title, fontsize=13, pad=10)
            ax.set_xlabel("x (lattice sites)")
            ax.set_ylabel("y (lattice sites)")
            ax.set_xticks([0, 20, 40]); ax.set_yticks([0, 20, 40])
            cb = fig.colorbar(image, ax=ax, shrink=.86, pad=.035)
            cb.ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
            cb.set_label("count per site", fontsize=11)
    fig.suptitle("CPM parcel-label trajectory", fontsize=19, y=.985)
    fig.text(.5, .94, "N = 40  |  42 initial parcels  |  seed 42  |  coupled MersenneTwister run",
             ha="center", fontsize=12)
    fig.text(.5, .118, "Maps: MCS 0-100, sampled every MCS; z = 20 (1-based). Grey: outside domain.",
             ha="center", fontsize=11)
    handles = [Patch(facecolor=COLOURS[i+1], label=name) for i, name in enumerate(SPECIES)]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .015), ncol=4,
               frameon=False, fontsize=10, columnspacing=1.2, handlelength=1.1)
    base = output / FIGURE
    fig.savefig(base.with_suffix(".pdf"), dpi=260)
    fig.savefig(base.with_suffix(".png"), dpi=220)
    plt.close(fig)
    extracted = subprocess.run(["pdftotext", "-layout", str(base.with_suffix(".pdf")), "-"],
                               check=True, capture_output=True, text=True).stdout
    if not extracted.strip() or "MCS 100" not in extracted:
        raise ValueError("Figure 5 text extraction is empty or lacks MCS 100")
    base.with_suffix(".txt").write_text(extracted)
    base.with_suffix(".sha256").write_text(sha(base.with_suffix(".pdf")) + "\n")
    manifest["figure5"] = {"name": FIGURE, "renderer": "Matplotlib Agg voxel-face rasterization",
                           "ray_tracing": False, "camera": {"elevation": 25, "azimuth": 35},
                           "central_plane_julia_z": 20,
                           "renderer_version": matplotlib.__version__,
                           "renderer_source_sha256": sha(Path(__file__)),
                           "source_analysis_sha256": sha(run / "analysis" / "label_dynamics.h5")}
    return base


def all_interior(summary, stride, kind):
    return next(row for row in summary["cadences"][str(stride)][kind] if row["stratum"] == "all_interior")


def results_tex(run: Path, manifest: dict) -> str:
    summary = json.loads((run / "analysis" / "label_dynamics.json").read_text())
    parcel = all_interior(summary, 1, "parcel")
    species = all_interior(summary, 1, "species")
    denominator = parcel["denominator_sites"]
    run_id = manifest["run_id"].replace("_", r"\_")
    event_returns = manifest["verification"]["replay"]["event_return_episodes"]
    lines = [r"% Generated by tools/render_label_trajectory.py from verified labelled outputs.",
             r"\subsection{Sampled Parcel-Label Trajectory}\label{sec:label_trajectory}",
             f"The separate diagnostic run \\texttt{{{run_id}}} uses the manuscript coupled",
             r"""configuration ($N=40$, six computational parcels per species, seed 42,
             \texttt{MersenneTwister}) and records 101 labelled states from MCS 0 through 100.""",
             f"Its initial registry and lattice both contain 42 parcels, with six for each of seven species.",
             f"The fixed interior contains {denominator:,} sites, including initially empty sites.",
             "At one-MCS sampling, the parcel-ID sequence contains "
             f"{parcel['transition_count']:,} adjacent-frame label transitions and "
             f"{parcel['return_episode_count']:,} completed return episodes; the species sequence contains "
             f"{species['transition_count']:,} transitions and {species['return_episode_count']:,} returns.",
             "A return is an observed departure from the MCS-0 reference label followed by a later arrival at that label.",
             f"Unchanged-at-every-saved-frame parcel identity holds at "
             f"{parcel['persistent_site_counts'][-1]:,}/{denominator:,} interior sites "
             f"({100*parcel['persistence_fraction'][-1]:.2f}\\%); this includes empty-site identity.",
             r"Initial-species and initial-parcel strata, their denominators and persistence curves are retained in the run data.",
             r"Binary species-indicator means report sampled occupancy; integer label codes are never averaged.",
             "", r"\begin{table}[htbp]", r"\centering\small",
             r"""\caption{Sampling the same trajectory at different MCS intervals. Counts are summed over the fixed interior.
             Hidden reversals count site--interval pairs whose coarse endpoint parcel IDs agree despite changes visible in the one-MCS sequence.}""",
             r"\label{tab:label_cadence}", r"\begin{tabular}{@{}rrrr@{}}\toprule",
             r"Interval (MCS) & Parcel transitions & Parcel returns & Hidden reversals \\\midrule"]
    for stride in (1, 2, 5, 10):
        row = all_interior(summary, stride, "parcel")
        lines.append(f"{stride} & {row['transition_count']:,} & {row['return_episode_count']:,} & {row['hidden_reversal_intervals']:,} " + r"\\")
    lines += [r"\bottomrule\end{tabular}\end{table}", "",
              "The ordered accepted-copy record contains "
              f"{manifest['accepted_copy_count']:,} events and replays the complete label-state hash at every MCS.",
              "It records " + f"{event_returns['parcel']:,} parcel-ID and {event_returns['species']:,} species return episodes at event resolution.",
              r"These event counts and the sampled counts have different temporal resolution; only the event record resolves ordering within a sweep.",
              r"The instrumented and uninstrumented complete states and RNG continuation agree exactly; one-MCS and 100-MCS windows also agree.",
              r"The original coupled runner agrees at MCS 0, 20, 40, 60, 80 and 100. A tenfold uptake-parameter perturbation leaves the labelled trajectory unchanged.",
              r"The radiodialysis basis acknowledgement is restricted to this label diagnostic: mobile and sorbed quantities remain blocked and are neither analysed nor displayed.",
              r"This is one deterministic computational trajectory, without an uncertainty estimate or biological-effect inference.",
              "", r"\begin{figure}[htbp]", r"\centering",
              r"\includegraphics[width=\linewidth]{fig5_label_trajectory}",
              r"""\caption{CPM parcel-label trajectory. Fixed-camera categorical species snapshots at MCS 0, 30 and 100
              (top) accompany central-plane parcel-ID transition and return counts sampled every MCS (bottom).
              The plane is the lower central slice, Julia $z=20$ on the $40^3$ grid; counts retain empty interior sites as label 0.
              Rendering uses rasterized voxel faces with ray tracing disabled. The evidence concerns computational parcel labels,
              not membranes, extracellular polymeric substances (EPS), or biological reconstitution.}""",
              r"\label{fig:label_trajectory}", r"\end{figure}", ""]
    return "\n".join(lines)


def stage_and_build(run: Path, manifest: dict, base: Path, *, install: bool) -> None:
    staging = run / "manuscript"
    if staging.exists():
        raise FileExistsError("refusing existing manuscript staging directory")
    staging.mkdir()
    shutil.copytree(ROOT / "preprint" / "figures", staging / "figures")
    for tex in (ROOT / "preprint").glob("*.tex"):
        shutil.copy2(tex, staging / tex.name)
    snippet = results_tex(run, manifest)
    (staging / "label_trajectory_results.tex").write_text(snippet)
    for suffix in (".pdf", ".png", ".txt", ".sha256"):
        shutil.copy2(base.with_suffix(suffix), staging / "figures" / (FIGURE + suffix))
    source = (staging / MANUSCRIPT).read_text()
    if r"\input{label_trajectory_results}" not in source:
        raise ValueError("manuscript must explicitly include the new results fragment")
    logpath = staging / "build.log"
    with logpath.open("w") as log:
        subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", MANUSCRIPT],
                       cwd=staging, stdout=log, stderr=subprocess.STDOUT, check=True)
    if logpath.stat().st_size == 0:
        raise ValueError("empty manuscript build log")
    latex_log = (staging / Path(MANUSCRIPT).with_suffix(".log")).read_text()
    if "There were undefined references" in latex_log or "multiply defined" in latex_log:
        raise ValueError("unresolved manuscript references")
    pdf = staging / Path(MANUSCRIPT).with_suffix(".pdf")
    if pdf.stat().st_size == 0:
        raise ValueError("empty compiled manuscript")
    text = subprocess.run(["pdftotext", "-layout", str(pdf), "-"], check=True, capture_output=True, text=True).stdout
    if "Figure 5:" not in text or "Sampled Parcel-Label Trajectory" not in text:
        raise ValueError("Figure 5 or generated results absent from compiled manuscript")
    pdf.with_suffix(".txt").write_text(text)
    manifest["manuscript"] = {"tex_sha256": sha(staging / MANUSCRIPT),
                               "results_tex_sha256": sha(staging / "label_trajectory_results.tex"),
                               "figure_pdf_sha256": sha(base.with_suffix(".pdf")),
                               "pdf_sha256": sha(pdf), "build_log_sha256": sha(logpath),
                               "visual_inspection": "pending"}
    if install:
        # Explicit --install copies only the new figure and its generated
        # fragment. Existing figure bytes are checked again after the copy.
        for suffix in (".pdf", ".png", ".txt", ".sha256"):
            shutil.copy2(base.with_suffix(suffix), ROOT / "preprint" / "figures" / (FIGURE + suffix))
        shutil.copy2(staging / "label_trajectory_results.tex", ROOT / "preprint" / "label_trajectory_results.tex")
    for rel, expected in manifest["figure1_hashes"].items():
        if sha(ROOT / rel) != expected:
            raise ValueError("Figure 1 changed")
    manifest["status"] = "manuscript_built"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--install", action="store_true", help="copy the new Figure 5 and results fragment into preprint")
    args = parser.parse_args()
    run = args.run_dir.resolve()
    manifest = verified_manifest(run)
    base = render(run, manifest)
    stage_and_build(run, manifest, base, install=args.install)
    save_manifest(run, manifest)
    print(run / "manuscript" / Path(MANUSCRIPT).with_suffix(".pdf"))


if __name__ == "__main__":
    main()
