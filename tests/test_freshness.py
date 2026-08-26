"""
Reload detection and freshness authority must not be the same signal.

The load-bearing test here is `test_touch_reloads_but_does_not_refresh`: a file
whose bytes changed and whose semantic timestamp did not must come back RELOADED
and still STALE. If mtime alone could promote a layer to nominal, `touch` would
be a way to lie to yourself about how current your data is.
"""

import json
import os
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from geolocator.freshness import (  # noqa: E402
    AUTHORED, FALLBACK, LIVE, MTIME_AUTHORITY, NOMINAL, REFERENCE, STALE,
    UNAVAILABLE, TrackedSource, parse_ts,
)


def iso(dt):
    return dt.isoformat()


def write(path, generated, items):
    path.write_text(json.dumps({"generated": generated, "items": items}))


def live_source(path):
    return TrackedSource(
        id="t", path=path, layer_class=LIVE,
        loader=lambda p: json.loads(p.read_text()),
        timestamp_of=lambda d: (d or {}).get("generated"),
        count_of=lambda d: len((d or {}).get("items", [])),
        degraded_after_s=90, stale_after_s=180,
        empty={"generated": None, "items": []},
    )


OLD = iso(datetime.now(timezone.utc) - timedelta(days=3))
NEW = iso(datetime.now(timezone.utc))


# ── the correction, stated as a test ─────────────────────────────────────────
def test_touch_reloads_but_does_not_refresh(tmp_path):
    f = tmp_path / "world.json"
    write(f, OLD, [1, 2, 3])
    src = live_source(f)
    before = src.state()
    assert before.status == STALE and before.count == 3

    time.sleep(0.01)
    # Bytes change, semantic timestamp does not -- the strongest form of the
    # test: it proves a reread happened AND that freshness did not move.
    write(f, OLD, [1, 2, 3, 4, 5])
    after = src.state()

    assert after.count == 5, "file was not rereread"
    assert after.fingerprint != before.fingerprint, "fingerprint did not move"
    assert after.reloaded_at != before.reloaded_at, "reload was not recorded"
    assert after.source_timestamp == before.source_timestamp
    assert after.status == STALE, "mtime must never promote a layer to nominal"


def test_bare_touch_changes_nothing_but_the_fingerprint(tmp_path):
    f = tmp_path / "world.json"
    write(f, OLD, [1, 2, 3])
    src = live_source(f)
    before = src.state()
    time.sleep(0.01)
    os.utime(f, None)  # a literal `touch`
    after = src.state()
    assert after.status == STALE == before.status
    assert after.source_timestamp == before.source_timestamp
    assert after.count == before.count


def test_a_genuinely_newer_source_timestamp_does_refresh(tmp_path):
    f = tmp_path / "world.json"
    write(f, OLD, [1])
    src = live_source(f)
    assert src.state().status == STALE
    write(f, NEW, [1])
    assert src.state().status == NOMINAL, "a real timestamp advance must count"


def test_atomic_replace_is_detected(tmp_path):
    """The lattice poller writes temp-then-rename; st_ino is what catches it."""
    f = tmp_path / "world.json"
    write(f, OLD, [1])
    src = live_source(f)
    assert src.state().count == 1
    tmp = tmp_path / "world.json.tmp"
    write(tmp, OLD, [1, 2, 3, 4])
    os.utime(tmp, (os.stat(f).st_atime, os.stat(f).st_mtime))  # same mtime
    os.replace(tmp, f)
    assert src.state().count == 4, "atomic replace slipped past the fingerprint"


# ── classes ──────────────────────────────────────────────────────────────────
def test_authored_has_no_authority_and_says_so(tmp_path):
    f = tmp_path / "hand.csv"
    f.write_text("name\na\nb\n")
    src = TrackedSource(id="hand", path=f, layer_class=AUTHORED,
                        loader=lambda p: p.read_text().splitlines())
    st = src.state()
    assert st.status == FALLBACK
    assert st.authority == MTIME_AUTHORITY
    assert "not a source timestamp" in st.note


def test_reference_shows_a_vintage_and_never_rots(tmp_path):
    f = tmp_path / "ref.json"
    write(f, iso(datetime.now(timezone.utc) - timedelta(days=400)), [1])
    src = TrackedSource(id="ref", path=f, layer_class=REFERENCE,
                        loader=lambda p: json.loads(p.read_text()),
                        timestamp_of=lambda d: d.get("generated"),
                        count_of=lambda d: len(d["items"]))
    st = src.state()
    assert st.status == NOMINAL, "a versioned dataset has a vintage, not an age"
    assert st.authority == "dataset-vintage" and st.vintage


def test_live_degrades_then_goes_stale(tmp_path):
    f = tmp_path / "w.json"
    for delta, expect in ((0, NOMINAL), (120, "degraded"), (600, STALE)):
        write(f, iso(datetime.now(timezone.utc) - timedelta(seconds=delta)), [1])
        assert live_source(f).state().status == expect, f"at {delta}s"


# ── failure modes ────────────────────────────────────────────────────────────
def test_missing_file_is_unavailable(tmp_path):
    src = live_source(tmp_path / "nope.json")
    st = src.state()
    assert st.status == UNAVAILABLE and st.count in (0, None)


def test_unparseable_file_is_unavailable_not_a_crash(tmp_path):
    f = tmp_path / "bad.json"
    f.write_text("{not json")
    st = live_source(f).state()
    assert st.status == UNAVAILABLE and "JSONDecodeError" in (st.error or "")


def test_parse_ts_accepts_date_only_and_offset_forms():
    assert parse_ts("2026-08-17").tzinfo is not None      # provenance.retrieved
    assert parse_ts("2026-08-26T21:25:09Z").tzinfo is not None
    assert parse_ts("2026-08-26T21:25:09+00:00").tzinfo is not None
    assert parse_ts(None) is None and parse_ts("") is None and parse_ts("junk") is None
