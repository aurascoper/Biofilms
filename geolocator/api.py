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
  worldgrid — gridded world energy model (joined from lattice board port 4199)
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
import json
import os
from pathlib import Path

from fastapi import FastAPI, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
DATA_CSV = REPO / "power_plant_database_global.csv"
SAMPLE_CSV = ROOT / "sample_power_plants.csv"
DATA_DIR = ROOT / "data"
# Energy world grid from the lattice board (port 4199)
ENERGY_WORLD_JSON = Path(
    os.environ.get(
        "ENERGY_WORLD_JSON",
        str(Path.home() / "Developer/live_trading/research/lattice_board/data/energy_world.json"),
    )
)

app = FastAPI(title="Primary Fuel Source Geolocator", version="2.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

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
}


def _f(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


# ── Power plants (WRI) ────────────────────────────────────────────────────────
def _load_plants() -> list[dict]:
    path = DATA_CSV if DATA_CSV.exists() else SAMPLE_CSV
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
    return plants


# ── Generic CSV layer loader (nuclear, battery, gridcoin, stars) ──────────────
def _load_layer_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
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
    return out


# ── Energy world grid (from lattice board port 4199) ─────────────────────────
def _load_worldgrid() -> list[dict]:
    """Convert energy_world.json 4°-bin grid into GeoJSON points."""
    if not ENERGY_WORLD_JSON.exists():
        return []
    try:
        d = json.loads(ENERGY_WORLD_JSON.read_text())
    except Exception:
        return []
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
        lon = -180.0 + (float(gx) + 0.5) * bin_deg
        lat = -90.0 + (float(gy) + 0.5) * bin_deg
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
    return out


# ── Layer registry ────────────────────────────────────────────────────────────
LAYER_FILES = {
    "nuclear": DATA_DIR / "nuclear_fuel_cycle.csv",
    "battery": DATA_DIR / "battery_fuel_cycle.csv",
    "gridcoin": DATA_DIR / "gridcoin_boinc.csv",
    "stars": DATA_DIR / "star_systems.csv",
}

PLANTS = _load_plants()
LAYERS: dict[str, list[dict]] = {
    "power": PLANTS,
    "nuclear": _load_layer_csv(LAYER_FILES["nuclear"]),
    "battery": _load_layer_csv(LAYER_FILES["battery"]),
    "gridcoin": _load_layer_csv(LAYER_FILES["gridcoin"]),
    "worldgrid": _load_worldgrid(),
    "stars": _load_layer_csv(LAYER_FILES["stars"]),
}
SOURCE = "WRI Global Power Plant Database" if DATA_CSV.exists() else "built-in sample"


def _to_geojson(items: list[dict]) -> dict:
    feats = []
    for p in items:
        color = LAYER_COLORS.get("power", {}).get(p["color_key"], "#BDC3C7")
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
    return {
        "status": "ok",
        "source": SOURCE,
        "plant_count": len(PLANTS),
        "layers": {k: len(v) for k, v in LAYERS.items()},
    }


@app.get("/api/layers")
def layers():
    return {
        "layers": [
            {
                "id": k,
                "count": len(v),
                "color_field": "primary_fuel" if k == "power" else "stage",
                "colors": LAYER_COLORS.get(k, {}),
            }
            for k, v in LAYERS.items()
        ]
    }


@app.get("/api/fuels")
def fuels(layer: str = Query("power")):
    items = LAYERS.get(layer, [])
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
    items = LAYERS.get(layer, [])
    total_cap = sum(p["capacity_mw"] for p in items)
    by_fuel: dict[str, float] = {}
    by_country: dict[str, int] = {}
    for p in items:
        by_fuel[p["color_key"]] = by_fuel.get(p["color_key"], 0.0) + p["capacity_mw"]
        by_country[p["country"]] = by_country.get(p["country"], 0) + 1
    return {
        "layer": layer,
        "source": SOURCE,
        "count": len(items),
        "total_capacity_mw": round(total_cap, 1),
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
    items = LAYERS.get(layer, [])
    out = []
    for p in items:
        if fuel and p["color_key"] != fuel:
            continue
        if p["capacity_mw"] < min_capacity:
            continue
        if max_capacity is not None and p["capacity_mw"] > max_capacity:
            continue
        if country and country.lower() not in p["country"].lower():
            continue
        out.append(p)
        if len(out) >= limit:
            break
    return _to_geojson(out)


@app.get("/")
def index():
    return FileResponse(ROOT / "static" / "index.html")
