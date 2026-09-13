"""Read-only research view. No lattice financial fields and no trading authority.

Legacy bars have no availability receipts and are ineligible. A separate,
prospective producer may call record_bar with an explicitly completed UTC bar;
the collector's observation time can never be backdated to the bar's date.
"""

import hashlib
import json
import math
import sqlite3
from datetime import datetime, timezone
from itertools import combinations, pairwise
from pathlib import Path

from geolocator.market import FUEL_TICKERS, MARKET_PROXY, MIN_OVERLAP, _pearson, _pvalue


def _utc(value: str) -> datetime:
    dt = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if dt.tzinfo is None or dt.utcoffset().total_seconds() != 0:
        raise ValueError("timestamp_must_be_utc")
    return dt


def _digest(value) -> str:
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def record_bar(
    path: Path,
    *,
    symbol: str,
    completed_at: str,
    close: str,
    source_url: str,
    now: datetime | None = None,
) -> None:
    """Prospective ingestion, deliberately separate from the GET endpoint.

    The declared completed_at must come from a verified provider convention.
    The current MEXC backfill does not call this function: its Day1 boundary
    convention has not been accepted for the new causal lane.
    """
    instant = now or datetime.now(timezone.utc)
    completed = _utc(completed_at)
    if instant.tzinfo is None or completed > instant:
        raise ValueError("bar_not_completed")
    if symbol not in {s for members in FUEL_TICKERS.values() for s in members} | {
        MARKET_PROXY
    }:
        raise ValueError("unknown_public_proxy")
    price = float(close)
    if not math.isfinite(price) or price <= 0 or not source_url.startswith("https://"):
        raise ValueError("invalid_bar_source")
    available = instant.astimezone(timezone.utc).isoformat()
    payload = {
        "symbol": symbol,
        "completedAt": completed.isoformat(),
        "close": close,
        "sourceUrl": source_url,
    }
    with sqlite3.connect(path) as con:
        con.execute(
            "CREATE TABLE IF NOT EXISTS research_bars (symbol TEXT, completed_at TEXT, available_at TEXT, close TEXT, source_url TEXT, digest TEXT, PRIMARY KEY(symbol,completed_at,available_at))"
        )
        con.execute(
            "INSERT INTO research_bars VALUES(?,?,?,?,?,?)",
            (
                symbol,
                completed.isoformat(),
                available,
                close,
                source_url,
                _digest(payload),
            ),
        )


def validate_family(document: dict) -> None:
    expected = set(combinations(sorted(FUEL_TICKERS), 2))
    actual = [(row["a"], row["b"]) for row in document["correlations"]]
    if (
        len(actual) != len(expected)
        or set(actual) != expected
        or document["familySize"] != len(expected)
    ):
        raise ValueError("incomplete_correlation_family")
    computed = sum(row["status"] == "available" for row in document["correlations"])
    if document["computedCount"] != computed:
        raise ValueError("correlation_count_mismatch")


def research_snapshot(path: Path, now: datetime | None = None) -> dict:
    instant = (now or datetime.now(timezone.utc)).astimezone(timezone.utc)
    pairs = list(combinations(sorted(FUEL_TICKERS), 2))
    document = {
        "schemaId": "biofilms.fuel-research.v1",
        "status": "unavailable",
        "reason": "availability_receipts_missing",
        "sourceType": "live",
        "sourceTimestamp": None,
        "availableAt": None,
        "fetchedAt": instant.isoformat(),
        "orderAuthority": False,
        "lineage": ["public-fuel-proxies"],
        "latticeEvidence": {
            "status": "unavailable",
            "reason": "strategy_lineage_not_declared",
        },
        "familySize": len(pairs),
        "computedCount": 0,
        "control": MARKET_PROXY,
        "baskets": {
            fuel: {"members": members, "returns": []}
            for fuel, members in FUEL_TICKERS.items()
        },
        "correlations": [
            {
                "a": a,
                "b": b,
                "status": "unavailable",
                "reason": "availability_receipts_missing",
                "n": 0,
            }
            for a, b in pairs
        ],
    }
    if not path.is_file():
        document["reason"] = "market_store_absent"
        return document
    con = None
    try:
        con = sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)
        con.row_factory = sqlite3.Row
        if not con.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='research_bars'"
        ).fetchone():
            return document
        raw = [
            dict(row)
            for row in con.execute(
                "SELECT * FROM research_bars ORDER BY symbol,completed_at,available_at"
            )
        ]
    except sqlite3.Error:
        document["reason"] = "market_store_unreadable"
        return document
    finally:
        if con is not None:
            con.close()
    closes = {}
    accepted = []
    try:
        for row in raw:
            completed, available = _utc(row["completed_at"]), _utc(row["available_at"])
            if not completed <= available <= instant:
                raise ValueError("future_or_noncausal_bar")
            price = float(row["close"])
            if not math.isfinite(price) or price <= 0:
                raise ValueError("invalid_bar_price")
            payload = {
                "symbol": row["symbol"],
                "completedAt": row["completed_at"],
                "close": row["close"],
                "sourceUrl": row["source_url"],
            }
            if _digest(payload) != row["digest"]:
                raise ValueError("bar_digest_mismatch")
            closes.setdefault(row["symbol"], {})[completed] = (price, available)
            accepted.append(row)
    except (ValueError, TypeError, KeyError) as exc:
        document["reason"] = str(exc)
        return document
    if not accepted:
        return document
    # Complete constituent coverage and consecutive daily bars. No silent
    # basket reweighting or multiday return labeled as one day's return.
    per_symbol = {}
    for symbol, rows in closes.items():
        per_symbol[symbol] = {}
        times = sorted(rows)
        for a, b in pairwise(times):
            if (b - a).total_seconds() == 86400:
                per_symbol[symbol][b] = (
                    math.log(rows[b][0] / rows[a][0]),
                    max(rows[a][1], rows[b][1]),
                )
    baskets = {}
    for fuel, members in {**FUEL_TICKERS, "__market__": [MARKET_PROXY]}.items():
        times = set.intersection(
            *(set(per_symbol.get(symbol, {})) for symbol in members)
        )
        baskets[fuel] = {
            t: (
                sum(per_symbol[s][t][0] for s in members) / len(members),
                max(per_symbol[s][t][1] for s in members),
            )
            for t in sorted(times)
        }
        if fuel != "__market__":
            document["baskets"][fuel]["returns"] = [
                {
                    "sourceTimestamp": t.isoformat(),
                    "availableAt": available.isoformat(),
                    "logReturn": str(value),
                    "complete": True,
                }
                for t, (value, available) in baskets[fuel].items()
            ]
    rows = []
    for a, b in pairs:
        shared = sorted(set(baskets[a]) & set(baskets[b]) & set(baskets["__market__"]))
        row = {
            "a": a,
            "b": b,
            "n": len(shared),
            "status": "unavailable",
            "reason": "insufficient_completed_overlap",
            "window": {"start": shared[0].isoformat(), "stop": shared[-1].isoformat()}
            if shared
            else None,
        }
        if len(shared) >= MIN_OVERLAP:
            x, y, z = (
                [baskets[fuel][t][0] for t in shared] for fuel in (a, b, "__market__")
            )
            xy, xz, yz = _pearson(x, y), _pearson(x, z), _pearson(y, z)
            if None not in (xy, xz, yz) and (1 - xz * xz) * (1 - yz * yz) > 0:
                partial = max(
                    -1.0,
                    min(1.0, (xy - xz * yz) / math.sqrt((1 - xz * xz) * (1 - yz * yz))),
                )
                row.update(
                    status="available",
                    reason=None,
                    rPartial=partial,
                    rRaw=xy,
                    p=_pvalue(partial, len(shared) - 1),
                )
            else:
                row["reason"] = "degenerate_series"
        rows.append(row)
    computed = sorted(
        (row for row in rows if row["status"] == "available"), key=lambda row: row["p"]
    )
    keep = max(
        (i for i, row in enumerate(computed, 1) if row["p"] <= 0.05 * i / len(pairs)),
        default=0,
    )
    for i, row in enumerate(computed, 1):
        row["significant"] = i <= keep
    document.update(
        correlations=rows,
        computedCount=len(computed),
        sourceTimestamp=max(row["completed_at"] for row in accepted),
        availableAt=max(row["available_at"] for row in accepted),
    )
    if all(baskets[fuel] for fuel in FUEL_TICKERS):
        age = (
            instant - min(max(baskets[fuel]) for fuel in FUEL_TICKERS)
        ).total_seconds()
        document.update(
            status="available" if 0 <= age <= 172800 else "unavailable",
            reason=None if 0 <= age <= 172800 else "completed_bars_stale",
        )
    document["provenance"] = {
        "schemaId": "neuralcompose.provenance-envelope.v1",
        "assertionKind": "derivedDeterministically",
        "method": {
            "methodId": "biofilms.fuel-research.v1",
            "softwareId": "biofilms-geolocator",
            "softwareVersion": "research-v1",
            "gitCommit": None,
            "parametersDigest": _digest(
                {
                    "baskets": FUEL_TICKERS,
                    "control": MARKET_PROXY,
                    "minimumOverlap": MIN_OVERLAP,
                    "bhQ": 0.05,
                }
            ),
        },
        "inputs": [
            {
                "resourceKind": "public-completed-bar-receipts",
                "sha256Hex": _digest(accepted),
                "locator": None,
            }
        ],
        "confidence": None,
        "comparisonEmbeddingSpace": None,
    }
    validate_family(document)
    return document
