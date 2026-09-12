#!/usr/bin/env python3
"""
geolocator/api.py — 3D geo-spatial energy-infrastructure locator API.

Replicates the "Primary Fuel Source Geolocator" Tableau workbook (Biofilms repo,
WRI Global Power Plant Database) as a live 3D globe API, expanded to multiple
energy layers:

  power     — WRI Global Power Plant Database (34,936 plants)
  nuclear   — nuclear fuel cycle (mines, enrichment, reprocessing, waste)
  battery   — battery fuel cycle (lithium/cobalt/nickel, gigafactories, recycling)
  gridcoin  — Gridcoin/BOINC energy grid (project servers)
  worldgrid — gridded world energy model (joined from the local lattice source)
  stars     — speculative star-system energy layer (exoplanet hosts)

Endpoints:
  GET /api/plants?layer=power&fuel=Coal&min_capacity=100&limit=5000
      -> GeoJSON FeatureCollection
  GET /api/layers
      -> available layers + counts + color keys
  GET /api/fuels?layer=power
      -> distinct fuel/stage values + counts + colors
  GET /api/stats?layer=power
      -> global aggregates
  GET /api/health
      -> liveness

Data: loads power_plant_database_global.csv if present (WRI/Kaggle schema),
else falls back to a built-in sample. Other layers load from geolocator/data/*.csv
and the energy world grid from the lattice board's energy_world.json.

Run:
  python3 -m uvicorn geolocator.api:app --host 127.0.0.1 --port 8000
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
from datetime import datetime, timezone
from pathlib import Path

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from geolocator.freshness import AUTHORED, LIVE, REFERENCE, TrackedSource
from geolocator import imagery
from geolocator.lattice import get_lattice

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
DATA_CSV = REPO / "power_plant_database_global.csv"
SAMPLE_CSV = ROOT / "sample_power_plants.csv"
DATA_DIR = ROOT / "data"
# Energy world grid from the local lattice source
ENERGY_WORLD_JSON = Path(
    os.environ.get(
        "ENERGY_WORLD_JSON",
        str(Path.home() / "Developer/live_trading/research/lattice_board/data/energy_world.json"),
    )
)
# Rev 6A physical-state overlay -- a symlink agri_yield_pipeline's export script retargets to
# its latest dated, content-hashed export each run. Read-only; never imported as code.
AGRI_OVERLAY_JSON = Path(
    os.environ.get(
        "AGRI_OVERLAY_JSON",
        str(Path.home() / "Developer/agri_yield_pipeline/data/overlay/latest.json"),
    )
)

app = FastAPI(title="Primary Fuel Source Geolocator", version="2.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


# ── caching policy, stated once ───────────────────────────────────────────────
#
#   /api/*   no-store    live and derived data. /api/lattice carries a redacted
#                        view of a private local source; it has no business on
#                        disk.
#                        The imagery raster is the one exception below.
#   /app/*   no-cache    versionless ES modules. May be STORED, must be
#                        REVALIDATED. Without this, StaticFiles sends only
#                        Last-Modified and ETag, and a browser may apply
#                        heuristic freshness and reuse a module blind --
#                        observed directly: after an edit, layerManager.js came
#                        back with transferSize 0 while its siblings took 304s,
#                        so the page ran a mixed old/new module graph and the
#                        symptom presented as an application bug.
#   /        no-cache    HTML that references those versionless modules.
#
# Nothing here is content-hashed, so nothing here earns a long immutable cache.
# The day an asset is served from a content-addressed URL, it wants
# `max-age=31536000, immutable` instead -- and that belongs in this table.
NO_STORE_PREFIX = "/api/"
REVALIDATE_PATHS = ("/app/", "/", "/legacy")
# Public imagery is large, unchanging and not sensitive; revalidating is enough,
# and a 304 on 209 KB is worth having.
REVALIDATE_EXCEPTIONS = ("/api/imagery/blue-marble",)


@app.middleware("http")
async def cache_policy(request, call_next):
    response = await call_next(request)
    path = request.url.path
    if path.startswith(REVALIDATE_EXCEPTIONS):
        response.headers["Cache-Control"] = "no-cache"
    elif path.startswith(NO_STORE_PREFIX):
        response.headers["Cache-Control"] = "no-store"
    elif path.startswith(REVALIDATE_PATHS) or path in REVALIDATE_PATHS:
        response.headers["Cache-Control"] = "no-cache"
    return response

# ── Fuel color palette (Tableau-style categorical ramp) ───────────────────────
FUEL_COLORS = {
    "Coal": "#8B4513",
    "Gas": "#E07B39",
    "Hydro": "#2E86AB",
    "Oil": "#C0392B",
    "Nuclear": "#7D3C98",
    "Solar": "#F4D03F",
    "Wind": "#5DADE2",
    "Biomass": "#27AE60",
    "Geothermal": "#E67E22",
    "Waste": "#95A5A6",
    "Petcoke": "#566573",
    "Wave": "#1ABC9C",
    "Storage": "#AF7AC5",
    "Cogeneration": "#D35400",
    "Wave and Tidal": "#1ABC9C",
    "Other": "#BDC3C7",
}

# Per-layer color maps (keyed by the color field for that layer)
LAYER_COLORS: dict[str, dict[str, str]] = {
    "power": FUEL_COLORS,
    "nuclear": {
        "uranium_mine": "#8B4513",
        "enrichment": "#7D3C98",
        "reprocessing": "#C0392B",
        "waste_storage": "#566573",
    },
    "battery": {
        "lithium": "#F4D03F",
        "cobalt": "#5DADE2",
        "nickel": "#27AE60",
        "gigafactory": "#E67E22",
        "recycling": "#1ABC9C",
    },
    "gridcoin": {
        "boinc_server": "#5D9BFF",
    },
    "worldgrid": FUEL_COLORS,
    "stars": {
        "dyson_swarm": "#FF6B6B",
        "fusion_core": "#FFD93D",
        "stellar_lens": "#6BCB77",
        "gridcoin_relay": "#4D96FF",
        "radio_silent": "#9AA7C0",
    },
    "agri_overlay": {
        # Same severity naming as agri_yield_pipeline's own field-stress seam
        # (src/stress_alerts.py:_severity -- "stress"/"warn"/"info").
        "Stress": "#C0392B",
        "Warn": "#E07B39",
        "Normal": "#27AE60",
        "No data": "#9AA7C0",
    },
}


def _f(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


# ── Power plants (WRI) ────────────────────────────────────────────────────────
def _load_plants(path: Path) -> dict:
    plants: list[dict] = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            lat = _f(row.get("latitude"))
            lon = _f(row.get("longitude"))
            if not (-90 <= lat <= 90) or not (-180 <= lon <= 180):
                continue
            plants.append(
                {
                    "name": row.get("name") or "",
                    "country": row.get("country") or row.get("country_long") or "",
                    "latitude": lat,
                    "longitude": lon,
                    "color_key": row.get("primary_fuel") or "Other",
                    "capacity_mw": _f(row.get("capacity_mw")),
                    "generation_gwh_2013": _f(row.get("generation_gwh_2013")),
                    "commissioning_year": row.get("commissioning_year") or "",
                    "owner": row.get("owner") or "",
                }
            )
    return {"items": plants, "meta": {}}


# ── Generic CSV layer loader (nuclear, battery, gridcoin, stars) ──────────────
def _load_layer_csv(path: Path) -> dict:
    out: list[dict] = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            lat = _f(row.get("latitude"))
            lon = _f(row.get("longitude"))
            if not (-90 <= lat <= 90) or not (-180 <= lon <= 180):
                continue
            out.append(
                {
                    "name": row.get("name") or "",
                    "country": row.get("country") or "",
                    "latitude": lat,
                    "longitude": lon,
                    "color_key": row.get("stage") or row.get("energy_class") or "Other",
                    "capacity_mw": _f(row.get("capacity_mw")),
                    "status": row.get("status") or "",
                    "source": row.get("source") or "",
                    "note": row.get("note") or "",
                    "extra": {
                        "ra": row.get("ra") or "",
                        "dec": row.get("dec") or "",
                        "distance_ly": row.get("distance_ly") or "",
                        "star_type": row.get("star_type") or "",
                    },
                }
            )
    return {"items": out, "meta": {}}


# ── Energy world grid (from the local lattice source) ────────────────────────
def _load_worldgrid(path: Path) -> dict:
    """Convert energy_world.json 4°-bin grid into GeoJSON points.

    `provenance` is carried out alongside the points because it holds this
    dataset's only semantic timestamp (`retrieved`). Dropping it here is what
    used to leave the filesystem mtime as the sole -- and misleading --
    freshness signal.
    """
    d = json.loads(path.read_text())
    nx = int(d.get("nx", 90))
    ny = int(d.get("ny", 45))
    bin_deg = float(d.get("bin_deg", 4.0))
    fuels = d.get("fuels", [])
    grid = d.get("grid", [])
    out: list[dict] = []
    for row in grid:
        if len(row) < 5:
            continue
        gx, gy, fuel_idx, mw, plants = row[0], row[1], row[2], row[3], row[4]
        # gy is stored north-to-south (row 0 = +90°), so latitude counts down.
        lon = -180.0 + (float(gx) + 0.5) * bin_deg
        lat = 90.0 - (float(gy) + 0.5) * bin_deg
        fuel = fuels[int(fuel_idx)] if 0 <= int(fuel_idx) < len(fuels) else "Other"
        out.append(
            {
                "name": f"grid cell ({gx},{gy})",
                "country": "",
                "latitude": lat,
                "longitude": lon,
                "color_key": fuel,
                "capacity_mw": _f(mw),
                "status": "operating",
                "source": "WRI energy world grid",
                "note": f"{int(plants)} plants in cell",
                "extra": {},
            }
        )
    return {"items": out, "meta": {"provenance": d.get("provenance") or {}}}


def _agri_stress_class(cell: dict) -> str:
    """Same info/warn/stress convention agri_yield_pipeline's own field-stress seam already
    uses (src/stress_alerts.py:_severity), applied to the county-level NDVI z here."""
    z = (cell.get("ndvi") or {}).get("z")
    if z is None:
        return "No data"
    az = abs(z)
    if az >= 2.0:
        return "Stress"
    if az >= 1.5:
        return "Warn"
    return "Normal"


def _load_agri_overlay(path: Path) -> dict:
    """Rev 6A agri_yield_pipeline physical-state overlay (see that repo's
    scripts/export_physical_overlay.py). An untrusted cross-repo file read -- hash-verifies the
    export before accepting it, same discipline as the corruption check in live_trading's
    energy_market_bridge_probe*.py snapshot loaders. Never imports agri_yield_pipeline's code,
    only reads its frozen JSON export.

    Verifies `payload_sha256` (schema_version>=2) over the whole payload minus that key itself --
    not just `cells`. A cells-only hash (schema_version 1's `cells_sha256`, kept here only for
    backward compatibility with exports written before this fix) would let a corrupted or
    hand-edited `generated_at`/`source_git_sha`/`provenance` pass untouched, since none of those
    live inside `cells`."""
    d = json.loads(path.read_text())
    cells = d.get("cells") or []

    payload_expected = d.get("payload_sha256")
    if payload_expected is not None:
        unhashed = {k: v for k, v in d.items() if k != "payload_sha256"}
        actual = hashlib.sha256(
            json.dumps(unhashed, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        if actual != payload_expected:
            raise ValueError(
                f"agri overlay corruption: payload hash mismatch "
                f"(recorded {payload_expected[:12]}..., actual {actual[:12]}...)"
            )
    else:
        # schema_version 1 fallback: only cells were ever hashed.
        cells_expected = d.get("cells_sha256")
        if cells_expected is not None:
            actual = hashlib.sha256(
                json.dumps(cells, sort_keys=True, separators=(",", ":")).encode("utf-8")
            ).hexdigest()
            if actual != cells_expected:
                raise ValueError(
                    f"agri overlay corruption: cells hash mismatch "
                    f"(recorded {cells_expected[:12]}..., actual {actual[:12]}...)"
                )

    items: list[dict] = []
    for cell in cells:
        lat, lon = cell.get("lat"), cell.get("lon")
        if lat is None or lon is None or not (-90 <= lat <= 90) or not (-180 <= lon <= 180):
            continue
        ndvi = cell.get("ndvi") or {}
        z = ndvi.get("z")
        items.append(
            {
                "name": cell.get("name") or "",
                "country": cell.get("state") or "",
                "latitude": lat,
                "longitude": lon,
                "color_key": _agri_stress_class(cell),
                # Reused as a numeric marker-size magnitude, not a true capacity -- this layer
                # has no MW figure. |NDVI z| keeps larger markers on the more anomalous
                # counties. None (not 0.0) when no NDVI baseline exists -- "unknown" and "zero
                # magnitude" are different claims, and /api/stats + /api/plants are the ones
                # responsible for handling an unknown capacity_mw correctly, not this loader
                # papering over it by picking a number.
                "capacity_mw": abs(z) if z is not None else None,
                "status": cell.get("kind") or "county",
                "source": "agri_yield_pipeline",
                "note": f"NDVI z={z:.2f}" if z is not None else "no NDVI baseline for this county",
                "extra": {
                    "ndvi": cell.get("ndvi"),
                    "ndvi_unavailable_reason": cell.get("ndvi_unavailable_reason"),
                    "weather": cell.get("weather"),
                    "yield_sensitivity": cell.get("yield_sensitivity"),
                    "sar_vv_db": cell.get("sar_vv_db"),
                },
            }
        )
    return {
        "items": items,
        "meta": {
            "provenance": d.get("provenance") or {},
            "source_git_sha": d.get("source_git_sha"),
            "generated_at": d.get("generated_at"),
            "county_set_note": d.get("county_set_note"),
        },
    }


def _agri_vintage(payload) -> str | None:
    return ((payload or {}).get("meta") or {}).get("generated_at")


# ── Layer registry ────────────────────────────────────────────────────────────
LAYER_FILES = {
    "nuclear": DATA_DIR / "nuclear_fuel_cycle.csv",
    "battery": DATA_DIR / "battery_fuel_cycle.csv",
    "gridcoin": DATA_DIR / "gridcoin_boinc.csv",
    "stars": DATA_DIR / "star_systems.csv",
}

POWER_CSV = DATA_CSV if DATA_CSV.exists() else SAMPLE_CSV


# The vendored WRI snapshot declares its own provenance beside it, so `power`
# no longer borrows a retrieval date from the derived grid and label it a
# vintage. The distinction is the whole point of the model: 1.3.0 is what the
# data depicts, 2026-08-17 is when the bytes arrived.
WRI_PROVENANCE_JSON = REPO / "power_plant_database_global.provenance.json"

WRI_PROVENANCE = TrackedSource(
    id="wri_provenance", path=WRI_PROVENANCE_JSON, layer_class=REFERENCE,
    loader=lambda p: json.loads(p.read_text()), empty={},
)


def _wri(field: str):
    try:
        return (WRI_PROVENANCE.payload() or {}).get(field)
    except Exception:
        return None


def _worldgrid_retrieved(payload) -> str | None:
    return ((payload or {}).get("meta", {}).get("provenance") or {}).get("retrieved")


# Reload is stat-driven (st_mtime_ns + st_size + st_ino). Freshness comes from
# inside the data. The four authored CSVs pass timestamp_of=None on purpose:
# they carry no semantic timestamp, so they report FALLBACK forever rather than
# letting an mtime impersonate one.
_EMPTY_LAYER = {"items": [], "meta": {}}
_layer_count = lambda d: len((d or {}).get("items", []))

SOURCES: dict[str, TrackedSource] = {
    "power": TrackedSource(
        id="power", path=POWER_CSV, layer_class=REFERENCE,
        loader=_load_plants, empty=_EMPTY_LAYER, count_of=_layer_count,
        vintage_of=lambda _p: _wri("vintage"),
        retrieved_of=lambda _p: _wri("retrieved"),
    ),
    "worldgrid": TrackedSource(
        # Recomputed from the same WRI release, so it carries the same vintage.
        id="worldgrid", path=ENERGY_WORLD_JSON, layer_class=REFERENCE,
        loader=_load_worldgrid, empty=_EMPTY_LAYER, count_of=_layer_count,
        vintage_of=lambda _p: _wri("vintage"),
        retrieved_of=_worldgrid_retrieved,
    ),
    "nuclear": TrackedSource(
        id="nuclear", path=LAYER_FILES["nuclear"], layer_class=AUTHORED,
        loader=_load_layer_csv, empty=_EMPTY_LAYER, count_of=_layer_count,
    ),
    "battery": TrackedSource(
        id="battery", path=LAYER_FILES["battery"], layer_class=AUTHORED,
        loader=_load_layer_csv, empty=_EMPTY_LAYER, count_of=_layer_count,
    ),
    "gridcoin": TrackedSource(
        id="gridcoin", path=LAYER_FILES["gridcoin"], layer_class=AUTHORED,
        loader=_load_layer_csv, empty=_EMPTY_LAYER, count_of=_layer_count,
    ),
    "stars": TrackedSource(
        id="stars", path=LAYER_FILES["stars"], layer_class=AUTHORED,
        loader=_load_layer_csv, empty=_EMPTY_LAYER, count_of=_layer_count,
    ),
    "agri_overlay": TrackedSource(
        # Rev 6A: county-level NDVI/weather/yield-sensitivity physical-state overlay from
        # agri_yield_pipeline. No trading claim; see that repo's plan for scope. REFERENCE, not
        # LIVE -- the vintage is the export's generated_at, never today's clock.
        id="agri_overlay", path=AGRI_OVERLAY_JSON, layer_class=REFERENCE,
        loader=_load_agri_overlay, empty=_EMPTY_LAYER, count_of=_layer_count,
        vintage_of=_agri_vintage,
        retrieved_of=_agri_vintage,
    ),
}
LAYER_IDS = tuple(SOURCES.keys())


def _max_bar_ts(path: Path) -> dict:
    """Newest daily bar in the cache. `ts` is epoch SECONDS."""
    import sqlite3
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        row = con.execute("SELECT MAX(ts), COUNT(*) FROM bars").fetchone()
    finally:
        con.close()
    if not row or row[0] is None:
        return {"max_ts": None, "rows": 0}
    return {"max_ts": int(row[0]), "rows": int(row[1])}


# The measured correlation band is only as current as its newest bar. Daily
# bars, so the thresholds are days: fine at one day old, degraded at two,
# stale past four. Nothing tells this layer it is behind -- it reports it.
# ── optional market adapter ───────────────────────────────────────────────────
#
# The measured correlation band needs a price source, and a price source is
# operator plumbing: venue API access plus a symbol universe. It is therefore an
# OPTIONAL ADAPTER rather than core code, and this deployment may legitimately
# ship without one.
#
# Contract — supply a module importable as `geolocator.market` exposing:
#
#     fuel_correlations() -> list[dict]
#         each dict: {a, b, r, r_raw, n, p, significant}
#         a, b            fuel names as they appear in FUEL_COLORS
#         r, r_raw        partial and raw correlation
#         n, p            sample size and p-value
#         significant     bool; only True entries are rendered
#
# Nothing else in this file needs to change to add one.
MARKET_ADAPTER_MODULE = "geolocator.market"


def market_adapter_installed() -> bool:
    try:
        return importlib.util.find_spec(MARKET_ADAPTER_MODULE) is not None
    except (ImportError, ValueError):
        return False


def market_correlation_state() -> dict:
    """State of the measured band's producer.

    An absent adapter must not read as "no correlations found". Absence of the
    producer and a null result from the producer are different facts, and
    collapsing them is the same error this codebase refuses everywhere else.
    """
    if not market_adapter_installed():
        return {
            "id": "market_correlation",
            "layer_class": LIVE,
            "status": "unavailable",
            "installed": False,
            "reason": f"optional local market adapter not installed ({MARKET_ADAPTER_MODULE})",
            "contract": "module exposing fuel_correlations() -> list[dict]",
        }
    state = MARKET_SOURCE.state().as_dict()
    state.update({"id": "market_correlation", "installed": True})
    if state.get("status") == "unavailable":
        state["reason"] = "market adapter installed, but its bar cache is empty or missing"
    return state


MARKET_SOURCE = TrackedSource(
    id="market_bars",
    path=Path(__file__).resolve().parent / "data" / "market_bars.sqlite",
    layer_class=LIVE,
    loader=_max_bar_ts,
    timestamp_of=lambda d: (
        datetime.fromtimestamp(d["max_ts"], tz=timezone.utc) if (d or {}).get("max_ts") else None
    ),
    degraded_after_s=2 * 86400,
    stale_after_s=4 * 86400,
    empty={"max_ts": None, "rows": 0},
    count_of=lambda d: (d or {}).get("rows", 0),
)
SOURCE = "WRI Global Power Plant Database" if DATA_CSV.exists() else "built-in sample"


def layer_items(layer: str) -> list[dict]:
    """Items for a KNOWN layer. Unknown ids 404 instead of returning [].

    Returning an empty FeatureCollection for a typo made "no such layer" and
    "no data in this layer" indistinguishable from the client.
    """
    src = SOURCES.get(layer)
    if src is None:
        raise HTTPException(
            status_code=404,
            detail=f"unknown layer {layer!r}; known layers: {', '.join(LAYER_IDS)}",
        )
    return (src.payload() or {}).get("items", [])


def _all_states() -> dict:
    return {k: v.state().as_dict() for k, v in SOURCES.items()}


def _to_geojson(items: list[dict], layer: str = "power") -> dict:
    """Palette is per layer. This took the "power" map for every layer, so every
    nuclear, battery, gridcoin and stars feature shipped grey while /api/fuels
    returned the correct colour, leaving the legend and the globe disagreeing.
    """
    palette = LAYER_COLORS.get(layer, {})
    feats = []
    for p in items:
        color = palette.get(p["color_key"], "#BDC3C7")
        feats.append(
            {
                "type": "Feature",
                "geometry": {
                    "type": "Point",
                    "coordinates": [p["longitude"], p["latitude"]],
                },
                "properties": {
                    "name": p["name"],
                    "country": p["country"],
                    "color_key": p["color_key"],
                    "capacity_mw": p["capacity_mw"],
                    "status": p.get("status", ""),
                    "source": p.get("source", ""),
                    "note": p.get("note", ""),
                    "extra": p.get("extra", {}),
                    "color": color,
                },
            }
        )
    return {"type": "FeatureCollection", "features": feats}


# ── Routes ────────────────────────────────────────────────────────────────────
@app.get("/api/health")
def health():
    states = _all_states()
    states["market_bars"] = MARKET_SOURCE.state().as_dict()
    states["market_correlation"] = market_correlation_state()
    states["blue_marble"] = imagery.state()
    return {
        "status": "ok",
        "source": SOURCE,
        "plant_count": len(layer_items("power")),
        "layers": {k: states[k].get("count") or 0 for k in LAYER_IDS},
        "freshness": states,
    }


@app.get("/api/layers")
def layers():
    states = _all_states()
    return {
        "layers": [
            {
                "id": k,
                "count": states[k].get("count") or 0,
                # worldgrid is coloured by fuel, not by stage -- it aliases FUEL_COLORS.
                # agri_overlay is coloured by NDVI stress class, not a fuel/stage/cycle
                # concept at all. A two-way ternary here has already once told the legend
                # something false for one layer; a dict keeps each layer's own answer.
                "color_field": {"power": "primary_fuel", "worldgrid": "primary_fuel",
                                 "agri_overlay": "stress_class"}.get(k, "stage"),
                "colors": LAYER_COLORS.get(k, {}),
                "layer_class": states[k].get("layer_class"),
                "status": states[k].get("status"),
                "authority": states[k].get("authority"),
                "source_timestamp": states[k].get("source_timestamp"),
                "vintage": states[k].get("vintage"),
                "retrieved_at": states[k].get("retrieved_at"),
                "note": states[k].get("note"),
            }
            for k in LAYER_IDS
        ]
    }


@app.get("/api/fuels")
def fuels(layer: str = Query("power")):
    items = layer_items(layer)
    counts: dict[str, int] = {}
    for p in items:
        counts[p["color_key"]] = counts.get(p["color_key"], 0) + 1
    colors = LAYER_COLORS.get(layer, {})
    return {
        "fuels": [
            {"fuel": k, "count": v, "color": colors.get(k, "#BDC3C7")}
            for k, v in sorted(counts.items(), key=lambda x: -x[1])
        ]
    }


@app.get("/api/stats")
def stats(layer: str = Query("power")):
    items = layer_items(layer)
    # capacity_mw is None for a legitimately-unknown value (e.g. agri_overlay counties with no
    # NDVI baseline) -- excluded from the sum rather than coerced to 0.0, which would make
    # "unknown" indistinguishable from "measured zero" in any total or average. The exclusion
    # count is always reported, not just implied by a total that's quietly smaller than count.
    known = [p for p in items if p["capacity_mw"] is not None]
    total_cap = sum(p["capacity_mw"] for p in known)
    by_fuel: dict[str, float] = {}
    by_country: dict[str, int] = {}
    for p in known:
        by_fuel[p["color_key"]] = by_fuel.get(p["color_key"], 0.0) + p["capacity_mw"]
    for p in items:
        by_country[p["country"]] = by_country.get(p["country"], 0) + 1
    return {
        "layer": layer,
        "source": SOURCE,
        "count": len(items),
        "total_capacity_mw": round(total_cap, 1),
        "unknown_capacity_count": len(items) - len(known),
        "capacity_by_fuel": {
            k: round(v, 1) for k, v in sorted(by_fuel.items(), key=lambda x: -x[1])
        },
        "top_countries": sorted(by_country.items(), key=lambda x: -x[1])[:10],
    }


@app.get("/api/plants")
def plants(
    layer: str = Query("power"),
    fuel: str | None = Query(None, description="Filter by color_key (fuel/stage)"),
    min_capacity: float = Query(0.0, ge=0.0),
    max_capacity: float | None = Query(None, ge=0.0),
    country: str | None = Query(None),
    limit: int = Query(5000, ge=1, le=50000),
):
    items = layer_items(layer)
    out = []
    excluded_unknown_capacity = 0
    for p in items:
        if fuel and p["color_key"] != fuel:
            continue
        cap = p["capacity_mw"]
        # An explicit min_capacity floor can't be verified against an unknown value, so it
        # excludes it -- conservative: don't claim something unproven clears a floor. The
        # *default* min_capacity=0.0 must not filter out unknowns just because 0.0 is a
        # technically-passable floor when cap happens to be known.
        if min_capacity > 0.0 and (cap is None or cap < min_capacity):
            if cap is None:
                excluded_unknown_capacity += 1
            continue
        if max_capacity is not None and (cap is None or cap > max_capacity):
            if cap is None:
                excluded_unknown_capacity += 1
            continue
        if country and country.lower() not in p["country"].lower():
            continue
        out.append(p)
        if len(out) >= limit:
            break
    geojson = _to_geojson(out, layer)
    # Always present, defaulting to 0 -- matching /api/stats's unknown_capacity_count, so a
    # consumer can tell "no filter active" apart from "filter active, nothing excluded" instead
    # of the two endpoints answering the same question in different shapes.
    geojson["excluded_unknown_capacity"] = excluded_unknown_capacity
    return geojson


# ── Correlation bands (links) ─────────────────────────────────────────────────
def _endpoint(pt: dict, frame: str) -> dict:
    """One end of a band, tagged with the frame it lives in.

    Terrestrial sites have lat/lon. Stars have no terrestrial position at all, so
    they carry ra/dec and the client resolves them onto the celestial shell. Fuels
    have no position of any kind, so they carry only a name.

    That last frame exists because a correlation between two price series has no
    location. Anchoring it at each fuel's largest cell put Gas and Oil 4 degrees
    apart in Japan; anchoring at capacity-weighted centroids put Gas on the
    Greenland ice sheet and Coal in Mongolia, where neither fuel exists. Both
    invent geography the statistic does not have, so the measured band gets none.
    """
    if frame == "fuel":
        return {"frame": "fuel", "name": pt}
    if frame == "celestial":
        extra = pt.get("extra") or {}
        return {
            "frame": "celestial",
            "name": pt["name"],
            "ra": _f(extra.get("ra")),
            "dec": _f(extra.get("dec")),
        }
    return {
        "frame": "geo",
        "name": pt["name"],
        "lat": pt["latitude"],
        "lon": pt["longitude"],
    }


def _band(a: dict, b: dict, band: str, status: str, color: str,
          label: str, stats: dict | None = None) -> dict:
    """A band. `status` is the epistemic tag and drives how the UI draws it.

    structural  - a real relationship in the data, no statistics claimed
    measured    - carries a statistic, and must ship the numbers behind it
    speculative - generated fiction, must never read as either of the above
    """
    both_geo = a["frame"] == "geo" and b["frame"] == "geo"
    props = {"band": band, "status": status, "color": color, "label": label,
             "a": a, "b": b}
    if stats:
        props["stats"] = stats
    return {
        "type": "Feature",
        "geometry": {
            "type": "LineString",
            "coordinates": [[a["lon"], a["lat"]], [b["lon"], b["lat"]]],
        } if both_geo else None,
        "properties": props,
    }


# Stage pairs that genuinely abut in the fuel cycle. All-pairs within a step is
# deliberate: uranium is a fungible global commodity, so "any mine can feed any
# enricher" is closer to true than naming one route. These are adjacency, not
# shipments, and the label says so.
#
# enrichment -> reprocessing is absent on purpose. Reprocessing takes irradiated
# fuel out of a reactor, and this layer has no reactor stage, so that edge would
# skip a step the data does not contain.
STRUCTURAL_STEPS = [
    ("nuclear", "uranium_mine", "enrichment", "#8B4513", "mined uranium feeds enrichment"),
    ("nuclear", "reprocessing", "waste_storage", "#C0392B", "reprocessing feeds waste storage"),
    ("battery", "lithium", "gigafactory", "#F4D03F", "lithium feeds cell manufacturing"),
    ("battery", "cobalt", "gigafactory", "#5DADE2", "cobalt feeds cell manufacturing"),
    ("battery", "nickel", "gigafactory", "#27AE60", "nickel feeds cell manufacturing"),
    ("battery", "gigafactory", "recycling", "#E67E22", "cells feed recycling"),
]


def _stage(layer: str, key: str) -> list[dict]:
    return [p for p in layer_items(layer) if p["color_key"] == key]


def _build_links() -> dict:
    bands: list[dict] = []

    # Structural: fuel-cycle stage adjacency.
    for layer, src, dst, color, label in STRUCTURAL_STEPS:
        for a in _stage(layer, src):
            for b in _stage(layer, dst):
                bands.append(_band(_endpoint(a, "geo"), _endpoint(b, "geo"),
                                   f"{layer}_chain", "structural", color, label))

    # Measured: fuel-basket return correlations, market-controlled. Carries no
    # geometry by design (see _endpoint). The client renders it as a paired
    # highlight of both fuels' cells, not as an arc. Only pairs surviving
    # Benjamini-Hochberg appear at all.
    try:
        from geolocator.market import fuel_correlations  # optional adapter
        for d in fuel_correlations():
            if not d.get("significant"):
                continue
            bands.append(_band(
                _endpoint(d["a"], "fuel"), _endpoint(d["b"], "fuel"),
                "fuel_correlation", "measured", "#5D9BFF",
                f"{d['a']} and {d['b']} returns move together",
                stats={"r_partial": d["r"], "r_raw": d["r_raw"], "n": d["n"],
                       "p": d["p"], "control": "QQQSTOCK_USDT",
                       "basis": "daily log returns, market-controlled"},
            ))
    except Exception:
        # Swallowed on purpose: a missing or failing adapter must not take the
        # structural and speculative bands with it. The ABSENCE is not silent --
        # market_correlation_state() reports it explicitly on /api/health and
        # /api/links, so "no adapter" never renders as "no correlations".
        pass

    # Speculative: the compute grid reaching toward the two nearest candidate
    # relay systems. Thematic, not derived from anything measured.
    relays = [p for p in layer_items("stars") if p["color_key"] == "gridcoin_relay"]
    for server in layer_items("gridcoin"):
        for star in relays:
            bands.append(_band(
                _endpoint(server, "geo"), _endpoint(star, "celestial"),
                "gridcoin_relay", "speculative", "#4D96FF",
                "SPECULATIVE: compute grid toward a candidate relay system",
            ))

    return {"type": "FeatureCollection", "features": bands}


_LINKS_CACHE: dict = {"key": None, "doc": None}

_LINK_SOURCES = ("nuclear", "battery", "gridcoin", "stars")


@app.get("/api/links")
def links():
    """Rebuilt only when a contributing source actually changed on disk."""
    key = tuple(SOURCES[k].state().fingerprint for k in _LINK_SOURCES)
    key += (MARKET_SOURCE.state().fingerprint,)
    if _LINKS_CACHE["key"] != key or _LINKS_CACHE["doc"] is None:
        _LINKS_CACHE["doc"] = _build_links()
        _LINKS_CACHE["key"] = key
    doc = dict(_LINKS_CACHE["doc"])
    doc["measured_freshness"] = market_correlation_state()
    return doc


@app.get("/api/lattice")
def lattice():
    """The only view of the upstream lattice this process ever serves.

    `view` is a server constant inside geolocator.lattice, never a parameter.
    The response is rebuilt from an allowlist rather than forwarded, so a field
    added upstream cannot become public data here by default.
    """
    return get_lattice()


@app.get("/api/imagery/blue-marble/meta")
def imagery_meta():
    """The imagery contract: provider, vintage, projection, attribution, and the
    SHA-256 of the exact bytes /api/imagery/blue-marble will serve."""
    return imagery.state()


@app.get("/api/imagery/blue-marble")
def imagery_raster():
    """Raster bytes only.

    There is deliberately no parameter here. The GetMap request lives in
    geolocator.imagery as a server constant, so no caller can steer the URL,
    the layer, the projection or the dimensions. An /api/imagery?url=... shape
    would be an open proxy.
    """
    st = imagery.state()
    path = imagery.raster_path()
    if path is None:
        raise HTTPException(status_code=503, detail=st.get("note") or "imagery unavailable")
    return FileResponse(
        path,
        media_type=imagery.MEDIA_TYPE,
        headers={
            "X-Imagery-Attribution": imagery.ATTRIBUTION,
            "X-Imagery-Vintage": imagery.DATASET_VINTAGE,
            "X-Imagery-Sha256": st.get("content_sha256") or "",
        },
    )


@app.get("/")
def index():
    return FileResponse(ROOT / "static" / "index.html")


@app.get("/legacy")
def legacy():
    """The pre-fusion globe, byte-for-byte, as a rollback path."""
    return FileResponse(ROOT / "static" / "index.legacy.html")


# The page is ES modules with no build step, so the browser fetches each one by
# path. Mounted LAST: a mount at /app would otherwise shadow nothing here, but
# keeping it after the routes makes the precedence obvious rather than lucky.
app.mount("/app", StaticFiles(directory=ROOT / "static" / "app"), name="app")
