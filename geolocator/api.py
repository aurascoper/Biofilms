#!/usr/bin/env python3
"""
geolocator/api.py — 3D geo-spatial power-plant locator API.

Replicates the "Primary Fuel Source Geolocator" Tableau workbook
(Biofilms repo, WRI Global Power Plant Database) as a live 3D globe API.

Endpoints:
  GET /api/plants?fuel=Coal&min_capacity=100&limit=5000
      -> GeoJSON FeatureCollection of power plants
  GET /api/fuels
      -> distinct primary_fuel values + counts
  GET /api/stats
      -> global aggregates (count, capacity by fuel, top countries)
  GET /api/health
      -> liveness

Data: loads power_plant_database_global.csv if present (WRI/Kaggle schema),
else falls back to a built-in sample so the globe renders immediately.

Run:
  uvicorn geolocator.api:app --host 127.0.0.1 --port 8000
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

from fastapi import FastAPI, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
DATA_CSV = REPO / "power_plant_database_global.csv"
SAMPLE_CSV = ROOT / "sample_power_plants.csv"

app = FastAPI(title="Primary Fuel Source Geolocator", version="1.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── Data loading ──────────────────────────────────────────────────────────────
# WRI Global Power Plant Database columns (Kaggle power_plant_database_global.csv)
WRI_FIELDS = {
    "name", "country", "country_long", "gppd_idnr", "capacity_mw",
    "latitude", "longitude", "primary_fuel", "other_fuel1", "other_fuel2",
    "other_fuel3", "commissioning_year", "owner", "source", "url",
    "geolocation_source", "wepp_id", "year_of_capacity_data",
    "generation_gwh_2013", "generation_gwh_2014", "generation_gwh_2015",
    "generation_gwh_2016", "generation_gwh_2017", "generation_gwh_2018",
    "generation_gwh_2019", "estimated_generation_gwh_2013",
    "estimated_generation_gwh_2014", "estimated_generation_gwh_2015",
    "estimated_generation_gwh_2016", "estimated_generation_gwh_2017",
    "estimated_generation_note_2013", "estimated_generation_note_2014",
    "estimated_generation_note_2015", "estimated_generation_note_2016",
    "estimated_generation_note_2017", "generation_data_source",
}

# Fuel color palette (matches a sensible Tableau-style categorical ramp)
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
    "Other": "#BDC3C7",
}


def _f(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


def _load_plants() -> list[dict]:
    """Load plants from the WRI CSV if present, else the built-in sample."""
    path = DATA_CSV if DATA_CSV.exists() else SAMPLE_CSV
    plants: list[dict] = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            lat = _f(row.get("latitude"))
            lon = _f(row.get("longitude"))
            if not (-90 <= lat <= 90) or not (-180 <= lon <= 180):
                continue
            cap = _f(row.get("capacity_mw"))
            gen = _f(row.get("generation_gwh_2013"))
            plants.append(
                {
                    "name": row.get("name") or "",
                    "country": row.get("country") or row.get("country_long") or "",
                    "latitude": lat,
                    "longitude": lon,
                    "primary_fuel": row.get("primary_fuel") or "Other",
                    "capacity_mw": cap,
                    "generation_gwh_2013": gen,
                    "commissioning_year": row.get("commissioning_year") or "",
                    "owner": row.get("owner") or "",
                }
            )
    return plants


PLANTS = _load_plants()
SOURCE = "WRI Global Power Plant Database" if DATA_CSV.exists() else "built-in sample"


def _to_geojson(plants: list[dict]) -> dict:
    feats = []
    for p in plants:
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
                    "primary_fuel": p["primary_fuel"],
                    "capacity_mw": p["capacity_mw"],
                    "generation_gwh_2013": p["generation_gwh_2013"],
                    "commissioning_year": p["commissioning_year"],
                    "owner": p["owner"],
                    "color": FUEL_COLORS.get(p["primary_fuel"], FUEL_COLORS["Other"]),
                },
            }
        )
    return {"type": "FeatureCollection", "features": feats}


# ── Routes ────────────────────────────────────────────────────────────────────
@app.get("/api/health")
def health():
    return {"status": "ok", "source": SOURCE, "plant_count": len(PLANTS)}


@app.get("/api/fuels")
def fuels():
    counts: dict[str, int] = {}
    for p in PLANTS:
        counts[p["primary_fuel"]] = counts.get(p["primary_fuel"], 0) + 1
    return {
        "fuels": [
            {"fuel": k, "count": v, "color": FUEL_COLORS.get(k, FUEL_COLORS["Other"])}
            for k, v in sorted(counts.items(), key=lambda x: -x[1])
        ]
    }


@app.get("/api/stats")
def stats():
    total_cap = sum(p["capacity_mw"] for p in PLANTS)
    by_fuel: dict[str, float] = {}
    by_country: dict[str, int] = {}
    for p in PLANTS:
        by_fuel[p["primary_fuel"]] = by_fuel.get(p["primary_fuel"], 0.0) + p["capacity_mw"]
        by_country[p["country"]] = by_country.get(p["country"], 0) + 1
    return {
        "source": SOURCE,
        "plant_count": len(PLANTS),
        "total_capacity_mw": round(total_cap, 1),
        "capacity_by_fuel": {
            k: round(v, 1) for k, v in sorted(by_fuel.items(), key=lambda x: -x[1])
        },
        "top_countries": sorted(by_country.items(), key=lambda x: -x[1])[:10],
    }


@app.get("/api/plants")
def plants(
    fuel: str | None = Query(None, description="Filter by primary_fuel"),
    min_capacity: float = Query(0.0, ge=0.0),
    max_capacity: float | None = Query(None, ge=0.0),
    country: str | None = Query(None),
    limit: int = Query(5000, ge=1, le=50000),
):
    out = []
    for p in PLANTS:
        if fuel and p["primary_fuel"] != fuel:
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
