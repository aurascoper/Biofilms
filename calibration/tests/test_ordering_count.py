"""The 256-seed ordering count that two figures print is the count in the committed ensemble.

Fig 2's label (biofilms_potts.jl) and panel A of the melanin ensemble figure (WIDER in
diagnostics/melanin_ensemble/figure.jl) both print how many seeds fail the ordering
C. sphaerospermum > C. neoformans > A. niger. Both were typed from ledger row SEED42-SCOPE-01
and nothing read them against an ensemble. The Julia model's own ensemble holds the answer:
radial_n40_p6_seeds42-297.csv, arm full, MCS 100, per-site mean_melanin.

Each figure is read twice, at its source string and at the .txt extracted from its PDF, so an
ensemble re-run that moves the count fails here until both figures are edited AND re-rendered.

WHAT BOUNDS IT. The failing count and the ensemble size, which is all either figure prints; not
which seeds fail, and not the range "42:297" printed beside the count in figure.jl. The
ordering is strict, so a tie fails it; the committed CSV has none.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
ENSEMBLE = Path("diagnostics/melanin_ensemble/radial_n40_p6_seeds42-297.csv")
ORDER = (3, 1, 5)   # species indices: C. sphaerospermum > C. neoformans > A. niger
PRINTED = (Path("biofilms_potts.jl"),
           Path("preprint/figures/fig2_melanin_accumulation.txt"),
           Path("diagnostics/melanin_ensemble/figure.jl"),
           Path("diagnostics/melanin_ensemble/melanin_ensemble_n40_p6.txt"))
COUNT = re.compile(r"fails in (\d+) of (\d+)")


def melanin_at_100(path: Path) -> dict[int, dict[int, float]]:
    """{seed: {species: mean_melanin}} for the full arm at MCS 100."""
    out: dict[int, dict[int, float]] = {}
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            if r["arm"] == "full" and r["mcs"] == "100":
                out.setdefault(int(r["seed"]), {})[int(r["species"])] = float(r["mean_melanin"])
    return out


def failing_seeds(path: Path) -> tuple[int, list[int]]:
    """The ensemble size, and the seeds whose ordering does not hold."""
    cs, cn, an = ORDER
    m = melanin_at_100(path)
    return len(m), sorted(s for s, v in m.items() if not v[cs] > v[cn] > v[an])


def problems(ensemble: Path, root: Path = REPO) -> list[tuple[Path, list, tuple]]:
    """Every printed location whose one count is not (failing, size). None found is a problem."""
    n, failing = failing_seeds(ensemble)
    want = (len(failing), n)
    out = []
    for rel in PRINTED:
        text = (root / rel).read_text(encoding="utf-8")
        found = [tuple(map(int, m.groups())) for m in COUNT.finditer(text)]
        if found != [want]:
            out.append((rel, found, want))
    return out


def test_both_figures_print_the_count_the_ensemble_gives():
    """THE ONE THAT MATTERS."""
    assert problems(REPO / ENSEMBLE) == []


def test_one_flipped_seed_is_reported_at_every_printed_location(tmp_path):
    source = (REPO / ENSEMBLE).read_text(encoding="utf-8").splitlines(keepends=True)
    n, failing = failing_seeds(REPO / ENSEMBLE)
    # The victim is taken from the file: the lowest seed whose ordering holds.
    victim = min(set(melanin_at_100(REPO / ENSEMBLE)) - set(failing))
    cs, cn, _ = ORDER

    def row(species):
        key = [str(victim), "full", "100", str(species)]
        hits = [i for i, line in enumerate(source) if line.split(",")[:4] == key]
        assert len(hits) == 1, (species, hits)
        return hits[0]

    i_cs, i_cn = row(cs), row(cn)
    mutated = list(source)
    head_cs, mel_cs = source[i_cs].rstrip("\n").rsplit(",", 1)
    head_cn, mel_cn = source[i_cn].rstrip("\n").rsplit(",", 1)
    mutated[i_cs] = f"{head_cs},{mel_cn}\n"
    mutated[i_cn] = f"{head_cn},{mel_cs}\n"
    scratch = tmp_path / ENSEMBLE.name
    scratch.write_text("".join(mutated), encoding="utf-8")

    # The mutation applied, and applied at the victim's two rows and nowhere else.
    assert [i for i in range(len(source)) if mutated[i] != source[i]] == sorted([i_cs, i_cn])
    assert failing_seeds(scratch) == (n, sorted(failing + [victim]))

    assert problems(scratch) == [(rel, [(len(failing), n)], (len(failing) + 1, n))
                                 for rel in PRINTED]


@pytest.mark.parametrize("victim", PRINTED, ids=[p.name for p in PRINTED])
def test_a_location_that_prints_no_count_is_refused(tmp_path, victim):
    for rel in PRINTED:
        (tmp_path / rel).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / rel).write_text((REPO / rel).read_text(encoding="utf-8"), encoding="utf-8")
    text = (tmp_path / victim).read_text(encoding="utf-8")
    [m] = list(COUNT.finditer(text))
    (tmp_path / victim).write_text(text[:m.start()] + text[m.end():], encoding="utf-8")
    assert COUNT.search((tmp_path / victim).read_text(encoding="utf-8")) is None

    n, failing = failing_seeds(REPO / ENSEMBLE)
    assert problems(REPO / ENSEMBLE, root=tmp_path) == [(victim, [], (len(failing), n))]
