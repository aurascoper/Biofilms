#!/usr/bin/env python3
"""
geolocator/lattice.py — the redaction boundary between the trading lattice and
this globe.

The lattice board on :4199 serves a live trading book. It offers a redacted
`?view=stream` mode, and this module asks for it -- but does not *trust* it.
`view=stream` is a property the lattice board establishes, and the lattice board
is free to change what it means without telling us. This process binds 0.0.0.0
with CORS "*", has no auth, and has already been exposed through a Cloudflare
tunnel, so a field added upstream tomorrow must not be able to become public
data by default.

So nothing is forwarded. A NEW document is built key by key from an allowlist:

  - no dict(upstream), no {**upstream}, no .update()
  - unknown keys are dropped by construction, not by filtering
  - values are coerced to scalars, so an allowlisted KEY cannot smuggle a
    nested object -- a key allowlist with unbounded values is not a boundary
  - `view` is a server constant and is never read from the client

What is deliberately NOT here: equity, available, margin_used, positions, and
by default has_position. The first four are already stripped upstream; they are
named here so that re-adding one is a visible edit to this list rather than an
accident.
"""

from __future__ import annotations

import json
import os
import urllib.error
import urllib.request

from geolocator.freshness import DEGRADED, NOMINAL, STALE, UNAVAILABLE, now_utc, parse_ts

# Same endpoint market.py already reads for its symbol universe. Read-only GET.
LATTICE_URL = os.environ.get("LATTICE_URL", "http://127.0.0.1:4199/api/lattice")

# A server constant. Not a parameter, not a default, not client-influenced.
UPSTREAM_VIEW = "stream"

TIMEOUT_S = float(os.environ.get("LATTICE_TIMEOUT_S", "3.0"))
CACHE_TTL_S = 2.0  # coalesces several browser tabs into one upstream GET

# ── the allowlist ─────────────────────────────────────────────────────────────
MISSION_FIELDS = ("phase", "phase_reason", "entries_allowed", "timestamp")

# `phase` and `phase_reason` are FREE TEXT written upstream ("growth",
# "manage-grow-above-floor"). A key allowlist cannot vet their contents: nothing
# stops a future reason string from reading "drawdown 87.5% from peak 9876".
# Since this module's whole premise is declining to trust upstream's judgment
# about what is safe, they are length-capped rather than passed through
# unbounded -- the same move serve.py makes when it truncates region-block
# messages to 160 chars. This BOUNDS the blast radius; it does not sanitise
# meaning, and it is not pretending to.
FREE_TEXT_FIELDS = ("phase", "phase_reason")
FREE_TEXT_CAP = 120
CELL_FIELDS = ("symbol", "quote", "class", "rank")
SIGNAL_FIELDS = ("direction", "score", "z_4h", "price", "liquidity", "asset")
ECON_GRIDCOIN_FIELDS = ("price_usd", "mcap_usd", "vol24h_usd", "chg24h_pct")

# has_position is NOT in CELL_FIELDS. Size is already hidden upstream, but
# presence still reveals which symbols are held, and this process is the wrong
# place for that to become a permanent public contract. Opt in explicitly.
EXPOSE_POSITION_PRESENCE = os.environ.get("GEOLOCATOR_EXPOSE_POSITION_PRESENCE") == "1"

# Live feed: the board is rewritten roughly every 60s by com.aurascoper.lattice-state.
DEGRADED_AFTER_S = 90.0
STALE_AFTER_S = 180.0

_cache: dict = {"at": None, "doc": None}


def _scalar(value):
    """Allowlisted keys carry scalars only. Anything else becomes None.

    Without this, `signal` could be allowlisted by name and still arrive as a
    nested object carrying whatever upstream decided to put there.
    """
    if isinstance(value, bool) or value is None:
        return value
    if isinstance(value, (int, float, str)):
        return value
    return None


def _cap(key, value):
    """Bound allowlisted free text. See FREE_TEXT_FIELDS."""
    if key in FREE_TEXT_FIELDS and isinstance(value, str) and len(value) > FREE_TEXT_CAP:
        return value[:FREE_TEXT_CAP] + "\u2026"
    return value


def _pick(src, fields) -> dict:
    """Build a new dict from `fields` only. Never mutates or copies `src`."""
    if not isinstance(src, dict):
        return {}
    return {k: _cap(k, _scalar(src.get(k))) for k in fields if k in src}


def _cell(raw) -> dict | None:
    if not isinstance(raw, dict):
        return None
    out = _pick(raw, CELL_FIELDS)
    if not out.get("symbol"):
        return None
    if EXPOSE_POSITION_PRESENCE:
        out["has_position"] = bool(raw.get("has_position"))
    sig = raw.get("signal")
    out["signal"] = _pick(sig, SIGNAL_FIELDS) if isinstance(sig, dict) else None
    return out


def _fetch_upstream() -> dict:
    url = f"{LATTICE_URL}?view={UPSTREAM_VIEW}"
    req = urllib.request.Request(url, headers={"User-Agent": "biofilms-geolocator/1.0"})
    with urllib.request.urlopen(req, timeout=TIMEOUT_S) as fh:
        return json.loads(fh.read().decode("utf-8"))


def _unavailable(reason: str) -> dict:
    return {
        "generated": None,
        "fetched_at": now_utc().isoformat(),
        "status": UNAVAILABLE,
        "authority": None,
        "age_s": None,
        "source": "lattice board :4199",
        "mission": None,
        "cells": [],
        "econ": None,
        "error": reason,
    }


def _build(upstream: dict) -> dict:
    state = upstream.get("state") if isinstance(upstream, dict) else None
    if not isinstance(state, dict):
        return _unavailable("upstream returned no state (poller may not have run)")

    generated = state.get("generated")
    semantic = parse_ts(generated)
    if semantic is None:
        status, age = UNAVAILABLE, None
    else:
        age = round((now_utc() - semantic).total_seconds(), 1)
        if age > STALE_AFTER_S:
            status = STALE
        elif age > DEGRADED_AFTER_S:
            status = DEGRADED
        else:
            status = NOMINAL

    raw_cells = state.get("cells")
    cells = []
    if isinstance(raw_cells, list):
        for raw in raw_cells:
            cell = _cell(raw)
            if cell is not None:
                cells.append(cell)

    econ_src = state.get("econ")
    grc = econ_src.get("gridcoin") if isinstance(econ_src, dict) else None
    econ = {"gridcoin": _pick(grc, ECON_GRIDCOIN_FIELDS)} if isinstance(grc, dict) else None

    return {
        "generated": _scalar(generated),
        "fetched_at": now_utc().isoformat(),
        "status": status,
        "authority": "source-timestamp" if semantic else None,
        "age_s": age,
        "source": "lattice board :4199",
        "mission": _pick(state.get("mission"), MISSION_FIELDS),
        "cells": cells,
        "econ": econ,
        "position_presence_exposed": EXPOSE_POSITION_PRESENCE,
    }


def get_lattice() -> dict:
    """The only thing the browser ever sees of the trading book."""
    now = now_utc()
    at = _cache["at"]
    if at is not None and (now - at).total_seconds() < CACHE_TTL_S:
        return _cache["doc"]
    try:
        doc = _build(_fetch_upstream())
    except urllib.error.URLError as exc:
        doc = _unavailable(f"lattice board unreachable: {exc.reason}")
    except (ValueError, TimeoutError, OSError) as exc:
        doc = _unavailable(f"{type(exc).__name__}: {exc}")
    _cache["at"] = now
    _cache["doc"] = doc
    return doc
