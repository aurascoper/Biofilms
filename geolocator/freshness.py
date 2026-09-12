#!/usr/bin/env python3
"""
geolocator/freshness.py — reload detection and freshness authority, kept apart.

Two questions that look alike and are not:

  WHEN TO REREAD                        IS IT CURRENT
  a stat fingerprint changed            a timestamp *inside the data* moved
  -> reload the file                    -> the data actually advanced

A filesystem mtime may trigger a reread. It may never, on its own, promote a
layer to nominal. `touch energy_world.json` reloads and leaves the status
exactly where it was, because touching a file does not make the world newer.

Every source therefore declares a CLASS, which decides where its freshness
authority comes from:

  live       expects continuous updates; measured against a clock.
             lattice -> state.generated, market bars -> MAX(bar_date).

  reference  a versioned dataset with a declared vintage. Never "stale by
             clock" -- 34,936 power plants retrieved on a date are not wrong
             a week later. Shows the vintage instead of an age.

  authored   hand-maintained, carries no semantic timestamp at all. Reports
             FALLBACK forever, labelled as a filesystem mtime so nobody reads
             it as a source timestamp.

That last case is the honest one: absence of an authority is a state to
report, not a gap to paper over with st_mtime.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

# The six feed states, ported from gods-eye-view's layerFeedState model.
NOMINAL = "nominal"
LOADING = "loading"
DEGRADED = "degraded"
STALE = "stale"
FALLBACK = "fallback"
UNAVAILABLE = "unavailable"

LIVE = "live"
REFERENCE = "reference"
AUTHORED = "authored"

MTIME_AUTHORITY = "filesystem-mtime"


def now_utc() -> datetime:
    return datetime.now(timezone.utc)


def parse_ts(value: Any) -> datetime | None:
    """Parse a semantic timestamp. Date-only values are accepted as UTC midnight.

    energy_world.json's provenance.retrieved is a bare date ("2026-08-17"); the
    lattice's state.generated is a full offset-aware ISO string. Both are valid
    authorities, so both parse here.
    """
    if not value:
        return None
    if isinstance(value, datetime):
        return value if value.tzinfo else value.replace(tzinfo=timezone.utc)
    text = str(value).strip()
    if not text:
        return None
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        dt = datetime.fromisoformat(text)
    except ValueError:
        return None
    return dt if dt.tzinfo else dt.replace(tzinfo=timezone.utc)


def _iso(dt: datetime | None) -> str | None:
    return dt.isoformat() if dt else None


@dataclass
class SourceState:
    """What the loader knows about one source, after a get()."""

    id: str
    layer_class: str
    status: str
    authority: str | None = None
    source_timestamp: str | None = None
    age_s: float | None = None
    vintage: str | None = None
    retrieved_at: str | None = None
    reloaded_at: str | None = None
    fingerprint: str | None = None
    count: int | None = None
    error: str | None = None
    note: str | None = None

    def as_dict(self) -> dict:
        return {k: v for k, v in self.__dict__.items() if v is not None}


@dataclass
class TrackedSource:
    """A file whose reload is stat-driven and whose freshness is data-driven.

    `loader`      Path -> payload.
    `timestamp_of` payload -> semantic timestamp, or None when the data carries
                  none. Returning None is a real answer, not a failure: it is
                  what puts a source into FALLBACK.
    """

    id: str
    path: Path
    layer_class: str
    loader: Callable[[Path], Any]
    timestamp_of: Callable[[Any], Any] | None = None
    count_of: Callable[[Any], int] | None = None
    # A reference dataset's authority is its DECLARED VINTAGE, which is often not
    # a date at all -- "1.3.0", "2004". Requiring a parseable timestamp here was
    # what pushed the WRI layer into borrowing its retrieval date and reporting
    # transport metadata under a `dataset-vintage` label.
    vintage_of: Callable[[Any], Any] | None = None
    retrieved_of: Callable[[Any], Any] | None = None
    vintage: str | None = None
    degraded_after_s: float = 90.0
    stale_after_s: float = 180.0
    empty: Any = field(default_factory=list)

    _payload: Any = field(default=None, init=False, repr=False)
    _fingerprint: tuple | None = field(default=None, init=False, repr=False)
    _loaded: bool = field(default=False, init=False, repr=False)
    _reloaded_at: datetime | None = field(default=None, init=False, repr=False)
    _error: str | None = field(default=None, init=False, repr=False)

    # ── reload detection ─────────────────────────────────────────────────────
    def _stat_fingerprint(self) -> tuple | None:
        """(st_mtime_ns, st_size, st_ino) -- an identity, not a clock.

        st_size and st_ino are in here because mtime alone misses an atomic
        replace that lands in the same nanosecond bucket, which is exactly how
        the lattice poller writes.
        """
        try:
            st = os.stat(self.path)
        except OSError:
            return None
        return (st.st_mtime_ns, st.st_size, st.st_ino)

    def _reload_if_changed(self) -> None:
        fp = self._stat_fingerprint()
        if fp is None:
            self._payload, self._fingerprint, self._loaded = self.empty, None, False
            self._error = f"missing: {self.path}"
            return
        if self._loaded and fp == self._fingerprint:
            return
        try:
            self._payload = self.loader(self.path)
            self._error = None
        except Exception as exc:  # a bad file must not take the process down
            self._payload = self.empty
            self._error = f"{type(exc).__name__}: {exc}"
        self._fingerprint = fp
        self._loaded = True
        self._reloaded_at = now_utc()

    # ── freshness authority ──────────────────────────────────────────────────
    def _resolve_state(self) -> SourceState:
        fp_repr = None if self._fingerprint is None else "%d/%d/%d" % self._fingerprint
        state = SourceState(
            id=self.id,
            layer_class=self.layer_class,
            status=LOADING,
            reloaded_at=_iso(self._reloaded_at),
            fingerprint=fp_repr,
            error=self._error,
        )
        try:
            state.count = self.count_of(self._payload) if self.count_of else len(self._payload)
        except (TypeError, AttributeError, KeyError):
            state.count = None

        if self._fingerprint is None:
            state.status = UNAVAILABLE
            return state
        if self._error:
            state.status = UNAVAILABLE
            return state

        if self.retrieved_of:
            got = self.retrieved_of(self._payload)
            state.retrieved_at = str(got) if got else None

        # A declared vintage settles a reference layer outright. No clock is
        # consulted, because a versioned snapshot does not drift toward wrong.
        if self.layer_class == REFERENCE:
            declared = self.vintage_of(self._payload) if self.vintage_of else self.vintage
            if declared:
                state.status = NOMINAL
                state.authority = "dataset-vintage"
                state.vintage = str(declared)
                return state

        semantic = parse_ts(self.timestamp_of(self._payload)) if self.timestamp_of else None

        if semantic is None:
            # No authority exists. Say so, and label the mtime for what it is.
            state.status = FALLBACK
            state.authority = MTIME_AUTHORITY
            state.note = "filesystem mtime - not a source timestamp"
            if self._fingerprint:
                mtime = datetime.fromtimestamp(self._fingerprint[0] / 1e9, tz=timezone.utc)
                state.source_timestamp = _iso(mtime)
            state.vintage = self.vintage
            return state

        state.source_timestamp = _iso(semantic)
        age = (now_utc() - semantic).total_seconds()
        state.age_s = round(age, 1)

        if self.layer_class == REFERENCE:
            # A versioned dataset has a vintage, not an age. It does not rot.
            state.status = NOMINAL
            state.authority = "dataset-vintage"
            state.vintage = self.vintage or _iso(semantic)
            return state

        state.authority = "source-timestamp"
        if age > self.stale_after_s:
            state.status = STALE
        elif age > self.degraded_after_s:
            state.status = DEGRADED
        else:
            state.status = NOMINAL
        return state

    # ── public ───────────────────────────────────────────────────────────────
    def get(self) -> tuple[Any, SourceState]:
        self._reload_if_changed()
        return self._payload, self._resolve_state()

    def payload(self) -> Any:
        return self.get()[0]

    def state(self) -> SourceState:
        return self.get()[1]
