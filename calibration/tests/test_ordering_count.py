"""The 256-seed ordering count that two figures print is the count in the committed ensemble.

Fig 2's label (biofilms_potts.jl) and panel A of the melanin ensemble figure (WIDER in
diagnostics/melanin_ensemble/figure.jl) both print how many seeds fail the ordering
C. sphaerospermum > C. neoformans > A. niger. Both were typed from ledger row SEED42-SCOPE-01
and nothing read them against an ensemble. The Julia model's own ensemble holds the answer:
radial_n40_p6_seeds42-297.csv, arm full, MCS 100, per-site mean_melanin.

Each figure is read twice, at its source string and at the .txt extracted from its PDF, so an
ensemble re-run that moves the count fails here until both figures are edited AND re-rendered.

WHAT BOUNDS IT. The failing count and the ensemble size, which is all either figure prints; not
which seeds fail. The CSV is held to the claim's own scope before it is counted: the seed set
is exactly 42..297, each seed carries exactly one full/MCS-100 row per producer species, and
each producer index carries the species_name the ordering claims. The literal "42:297" printed
beside the count in figure.jl is still read by nothing; what is checked is the CSV's seed set.
The ordering is strict, so a tie fails it; the committed CSV has none.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
ENSEMBLE = Path("diagnostics/melanin_ensemble/radial_n40_p6_seeds42-297.csv")
ORDER = (3, 1, 5)   # species indices: C. sphaerospermum > C. neoformans > A. niger
# What each ORDER index must be named in the CSV's species_name column, which radial_sweep.jl
# writes from SPECIES_NAMES in biofilms_potts.jl. An index remapped there changes the count's
# meaning without changing its arithmetic, so the loader refuses a row whose name disagrees.
NAMES = {3: "C. sphaerospermum", 1: "C. neoformans", 5: "A. niger"}
# The claim both figures print is scoped "Over seeds 42:297" and the CSV is named seeds42-297,
# so the expected seed set is that range, not whatever seeds the file happens to hold: 256
# distinct seeds with one of them swapped would otherwise count as the same ensemble.
SEEDS = range(42, 298)
PRINTED = (Path("biofilms_potts.jl"),
           Path("preprint/figures/fig2_melanin_accumulation.txt"),
           Path("diagnostics/melanin_ensemble/figure.jl"),
           Path("diagnostics/melanin_ensemble/melanin_ensemble_n40_p6.txt"))
COUNT = re.compile(r"fails in (\d+) of (\d+)")


def melanin_at_100(path: Path) -> dict[int, dict[int, float]]:
    """{seed: {species: mean_melanin}} for the full arm at MCS 100.

    Refuses a file that is not the ensemble the claim is about: a producer row whose
    species_name is not the name ORDER stands for, a duplicate (seed, species) row, a seed
    set other than exactly 42..297, or a seed missing a producer row.
    """
    out: dict[int, dict[int, float]] = {}
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            if r["arm"] != "full" or r["mcs"] != "100":
                continue
            seed, species = int(r["seed"]), int(r["species"])
            if species in NAMES and r["species_name"] != NAMES[species]:
                raise ValueError(f"species {species} is named {r['species_name']!r},"
                                 f" not {NAMES[species]!r}")
            if species in out.setdefault(seed, {}):
                raise ValueError(f"duplicate row for seed {seed} species {species}")
            out[seed][species] = float(r["mean_melanin"])
    unexpected, missing = sorted(set(out) - set(SEEDS)), sorted(set(SEEDS) - set(out))
    if unexpected or missing:
        raise ValueError(f"seed set is not {SEEDS}:"
                         f" unexpected {unexpected}, missing {missing}")
    short = sorted((s, sp) for s, v in out.items() for sp in ORDER if sp not in v)
    if short:
        raise ValueError(f"seeds missing a producer row: {short}")
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


# Negative controls for the loader's refusals. Each builds its known-bad from the committed
# CSV, shows a faithful copy still gives the committed answer, shows its mutation applied and
# where, and then requires the refusal that names the offender.


def _source() -> list[str]:
    return (REPO / ENSEMBLE).read_text(encoding="utf-8").splitlines(keepends=True)


def _rows(source: list[str], seed: int, species: int | None = None) -> list[int]:
    """Line indices of a seed's full/MCS-100 rows, optionally of one species."""
    def match(line: str) -> bool:
        f = line.split(",")
        return f[:3] == [str(seed), "full", "100"] and (species is None or f[3] == str(species))
    return [i for i, line in enumerate(source) if match(line)]


def _faithful_passes(tmp_path, source: list[str]) -> None:
    faithful = tmp_path / "faithful.csv"
    faithful.write_text("".join(source), encoding="utf-8")
    assert failing_seeds(faithful) == failing_seeds(REPO / ENSEMBLE)


def test_a_renumbered_seed_keeping_256_distinct_seeds_is_refused(tmp_path):
    """The completeness gap: a wrong seed set whose distinct-seed count is still 256."""
    source = _source()
    _faithful_passes(tmp_path, source)
    victim = min(melanin_at_100(REPO / ENSEMBLE))
    intruder = max(SEEDS) + 1
    rows = _rows(source, victim)
    assert rows, victim
    mutated = list(source)
    for i in rows:
        mutated[i] = f"{intruder},{source[i].split(',', 1)[1]}"
    assert [i for i in range(len(source)) if mutated[i] != source[i]] == rows
    scratch = tmp_path / ENSEMBLE.name
    scratch.write_text("".join(mutated), encoding="utf-8")
    with pytest.raises(ValueError,
                       match=re.escape(f"unexpected [{intruder}], missing [{victim}]")):
        failing_seeds(scratch)


def test_a_duplicated_producer_row_is_refused(tmp_path):
    source = _source()
    _faithful_passes(tmp_path, source)
    victim = min(melanin_at_100(REPO / ENSEMBLE))
    [i] = _rows(source, victim, ORDER[0])
    mutated = source[:i + 1] + [source[i]] + source[i + 1:]
    assert len(mutated) == len(source) + 1 and mutated[i] == mutated[i + 1] == source[i]
    scratch = tmp_path / ENSEMBLE.name
    scratch.write_text("".join(mutated), encoding="utf-8")
    with pytest.raises(ValueError,
                       match=re.escape(f"duplicate row for seed {victim} species {ORDER[0]}")):
        failing_seeds(scratch)


def test_a_deleted_producer_row_is_refused(tmp_path):
    """Named refusal, not the bare KeyError the count would otherwise raise."""
    source = _source()
    _faithful_passes(tmp_path, source)
    victim = min(melanin_at_100(REPO / ENSEMBLE))
    [i] = _rows(source, victim, ORDER[0])
    mutated = source[:i] + source[i + 1:]
    assert len(mutated) == len(source) - 1 and source[i] not in mutated
    scratch = tmp_path / ENSEMBLE.name
    scratch.write_text("".join(mutated), encoding="utf-8")
    with pytest.raises(ValueError,
                       match=re.escape(f"seeds missing a producer row: [({victim}, {ORDER[0]})]")):
        failing_seeds(scratch)


def test_a_remapped_species_name_is_refused(tmp_path):
    """The mapping gap: an index whose species_name is no longer what ORDER claims."""
    source = _source()
    _faithful_passes(tmp_path, source)
    victim = min(melanin_at_100(REPO / ENSEMBLE))
    [i] = _rows(source, victim, ORDER[0])
    fields = source[i].split(",")
    assert fields[4] == NAMES[ORDER[0]]   # the victim field holds the name being remapped
    fields[4] = NAMES[ORDER[1]]
    mutated = list(source)
    mutated[i] = ",".join(fields)
    assert [j for j in range(len(source)) if mutated[j] != source[j]] == [i]
    scratch = tmp_path / ENSEMBLE.name
    scratch.write_text("".join(mutated), encoding="utf-8")
    with pytest.raises(ValueError, match=re.escape(
            f"species {ORDER[0]} is named {NAMES[ORDER[1]]!r}, not {NAMES[ORDER[0]]!r}")):
        failing_seeds(scratch)
