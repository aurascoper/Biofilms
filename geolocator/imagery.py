#!/usr/bin/env python3
"""
geolocator/imagery.py — NASA GIBS base imagery, as a reference layer.

WHY THIS IS A LAYER AND NOT A GLOBE MATERIAL
Imagery arrives with a provider, a dataset vintage, a licence and an attribution
requirement, exactly like every other source here. Special-casing it into the
sphere's material would put the one dataset with a legal credit obligation
outside the model that tracks credit obligations. So it registers like the rest
and reports the same shape.

THE SEPARATION THAT MATTERS
    datasetVintage = 2004    epistemic. What the pixels depict.
    retrievedAt    = today   transport. When the bytes arrived.

Neither is staleness. Blue Marble Next Generation is a 2004 monthly composite;
it is twenty-two years old and completely current as a reference basemap,
because a reference dataset has a vintage rather than an age. STALE stays
reserved for live/measured sources whose age changes their validity -- which is
why market_bars can be STALE at 7.9 days while this is NOMINAL at 22 years.

If GIBS is unreachable and no cache exists, that is UNAVAILABLE, never STALE.
Not knowing and being out of date are different failures.

SSRF
The GetMap request is a server constant. There is no parameter, on any route,
that lets a caller influence the URL, the layer, the projection or the
dimensions -- the same rule that makes `view=stream` server-owned for the
lattice. An /api/imagery?url=... shaped endpoint would be an open proxy.
"""

from __future__ import annotations

import hashlib
import json
import os
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path

from geolocator.freshness import NOMINAL, REFERENCE, UNAVAILABLE, now_utc

CACHE_DIR = Path(
    os.environ.get("GEOLOCATOR_IMAGERY_CACHE",
                   str(Path(__file__).resolve().parent / "cache" / "imagery"))
)

# ── the request, frozen ───────────────────────────────────────────────────────
GIBS_ENDPOINT = "https://gibs.earthdata.nasa.gov/wms/epsg4326/best/wms.cgi"

# WMS 1.3.0 with EPSG:4326 is lat/lon axis order, so BBOX reads min-lat, min-lon,
# max-lat, max-lon. Getting this backwards yields a valid-looking sideways Earth.
GIBS_PARAMS = {
    "SERVICE": "WMS",
    "VERSION": "1.3.0",
    "REQUEST": "GetMap",
    "LAYERS": "BlueMarble_NextGeneration",
    "CRS": "EPSG:4326",
    "BBOX": "-90,-180,90,180",
    "WIDTH": "2048",
    "HEIGHT": "1024",
    "FORMAT": "image/jpeg",
    "STYLES": "",
}

DATASET_VINTAGE = "2004"
ATTRIBUTION = "NASA Earth Observatory"
MEDIA_TYPE = "image/jpeg"

TIMEOUT_S = float(os.environ.get("GEOLOCATOR_IMAGERY_TIMEOUT_S", "20"))
# A 2004 composite does not change. The TTL exists to re-verify the transport,
# not to chase updates.
CACHE_TTL = timedelta(days=30)
MAX_BYTES = 12 * 1024 * 1024

RASTER = CACHE_DIR / "blue_marble.jpg"
META = CACHE_DIR / "blue_marble.meta.json"


def request_descriptor() -> str:
    """The exact request, recorded so provenance is reproducible."""
    return f"{GIBS_ENDPOINT}?{urllib.parse.urlencode(GIBS_PARAMS)}"


def _fetch() -> bytes:
    req = urllib.request.Request(
        request_descriptor(), headers={"User-Agent": "biofilms-geolocator/1.0"})
    with urllib.request.urlopen(req, timeout=TIMEOUT_S) as fh:
        ctype = (fh.headers.get("Content-Type") or "").split(";")[0].strip()
        data = fh.read(MAX_BYTES + 1)
    if len(data) > MAX_BYTES:
        raise ValueError(f"imagery exceeds {MAX_BYTES} bytes")
    if ctype != MEDIA_TYPE:
        # GIBS reports errors as XML with a 200, so content-type is the real check.
        raise ValueError(f"unexpected content-type {ctype!r} (wanted {MEDIA_TYPE})")
    if not data.startswith(b"\xff\xd8\xff"):
        raise ValueError("response is not a JPEG")
    return data


def _read_meta() -> dict | None:
    try:
        return json.loads(META.read_text())
    except (OSError, ValueError):
        return None


def _write_cache(data: bytes) -> dict:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    meta = {
        "provider": "NASA GIBS",
        "dataset": GIBS_PARAMS["LAYERS"],
        "dataset_vintage": DATASET_VINTAGE,
        "projection": GIBS_PARAMS["CRS"],
        "service": GIBS_PARAMS["SERVICE"],
        "attribution": ATTRIBUTION,
        "media_type": MEDIA_TYPE,
        "width": int(GIBS_PARAMS["WIDTH"]),
        "height": int(GIBS_PARAMS["HEIGHT"]),
        "bytes": len(data),
        "content_sha256": hashlib.sha256(data).hexdigest(),
        "retrieved_at": now_utc().isoformat(),
        "request_descriptor": request_descriptor(),
    }
    # Write raster first: a meta file that describes bytes not yet on disk is
    # worse than a raster with no meta, which simply reads as a cache miss.
    tmp = RASTER.with_suffix(".jpg.tmp")
    tmp.write_bytes(data)
    os.replace(tmp, RASTER)
    META.write_text(json.dumps(meta, indent=2))
    return meta


def _cache_age() -> timedelta | None:
    meta = _read_meta()
    if not meta or not RASTER.exists():
        return None
    try:
        got = datetime.fromisoformat(meta["retrieved_at"])
    except (KeyError, ValueError):
        return None
    if not got.tzinfo:
        got = got.replace(tzinfo=timezone.utc)
    return now_utc() - got


def ensure() -> dict:
    """Resolve imagery, preferring cache within TTL, then GIBS, then stale cache.

    Never raises. The globe must boot whether or not NASA is up.
    """
    age = _cache_age()
    if age is not None and age < CACHE_TTL:
        return {"availability": "cached", "meta": _read_meta(),
                "note": "served from cache within TTL", "error": None}
    try:
        meta = _write_cache(_fetch())
        return {"availability": "live", "meta": meta, "note": None, "error": None}
    except (urllib.error.URLError, ValueError, TimeoutError, OSError) as exc:
        err = f"{type(exc).__name__}: {getattr(exc, 'reason', exc)}"
        if age is not None:
            return {"availability": "cached", "meta": _read_meta(),
                    "note": "provider unreachable; serving cached raster", "error": err}
        return {"availability": "unavailable", "meta": None,
                "note": "no imagery: provider unreachable and no cache", "error": err}


def raster_path() -> Path | None:
    return RASTER if RASTER.exists() else None


def state() -> dict:
    """Layer state in the same shape every other source reports.

    `source_timestamp` is None on purpose. This dataset has a vintage, not an
    observation time, and `retrieved_at` is transport metadata -- surfacing it
    as the source timestamp is exactly the confusion this model exists to avoid.
    """
    res = ensure()
    meta = res["meta"] or {}
    available = res["availability"] in ("live", "cached")
    return {
        "id": "blue_marble",
        "layer_class": REFERENCE,
        "kind": "imagery",
        "status": NOMINAL if available else UNAVAILABLE,
        "authority": "dataset-vintage" if available else None,
        "source_timestamp": None,
        "vintage": meta.get("dataset_vintage") if available else None,
        "age_s": None,
        "availability": res["availability"],
        "attribution": ATTRIBUTION,
        "provider": meta.get("provider"),
        "dataset": meta.get("dataset"),
        "projection": meta.get("projection"),
        "width": meta.get("width"),
        "height": meta.get("height"),
        "content_sha256": meta.get("content_sha256"),
        "retrieved_at": meta.get("retrieved_at"),
        "request_descriptor": meta.get("request_descriptor"),
        "note": res["note"],
        "error": res["error"],
    }
