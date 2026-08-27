#!/usr/bin/env python3
"""
geolocator/market.py — daily market bars for the correlation bands.

The globe's own generation series (WRI generation_gwh_2013..2019) has seven
annual points. Across thousands of candidate pairs that is not enough degrees of
freedom for any single band to survive multiple-testing correction, so the
correlation the globe draws is sourced here instead, where bars run to hundreds
or thousands of samples.

Read-only with respect to the trading repo: the symbol universe is pulled from
the lattice board over HTTP and the bars come from MEXC's public futures kline
endpoint. Nothing is written outside this repo.

Backfill (idempotent, safe to re-run):
  python3 -m geolocator.market backfill
"""

from __future__ import annotations

import json
import sqlite3
import sys
import time
import urllib.request
from pathlib import Path

DB_PATH = Path(__file__).resolve().parent / "data" / "market_bars.sqlite"
LATTICE_URL = "http://127.0.0.1:4199/api/lattice"
KLINE_URL = "https://futures.mexc.com/api/v1/contract/kline/{symbol}?interval=Day1"

SCHEMA = """
CREATE TABLE IF NOT EXISTS bars (
    symbol    TEXT    NOT NULL,
    ts        INTEGER NOT NULL,
    close     REAL    NOT NULL,
    quote_vol REAL,
    PRIMARY KEY (symbol, ts)
);
"""


def connect() -> sqlite3.Connection:
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(DB_PATH)
    con.executescript(SCHEMA)
    return con


def _get_json(url: str, timeout: int = 20) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "biofilms-geolocator/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as fh:
        return json.loads(fh.read().decode("utf-8"))


def universe() -> list[str]:
    """Symbol list from the lattice board. Read-only GET."""
    d = _get_json(LATTICE_URL)
    cells = ((d.get("state") or {}).get("cells")) or []
    return [c["symbol"] for c in cells if c.get("symbol")]


def fetch_bars(symbol: str) -> list[tuple[int, float, float]]:
    """Daily bars for one perp -> [(ts, close, quote_vol)].

    The futures kline endpoint returns parallel arrays (time[], close[], vol[]),
    not array-of-arrays like spot, and needs no start param: it hands back the
    full listed history, capped near 2000 bars. vol is base-coin, so quote volume
    is close * vol.

    Day1 is deliberate. Most of this universe is equity-linked perps that track US
    market hours while crypto trades around the clock, so intraday bars would fold
    session boundaries into any correlation drawn from them.
    """
    d = _get_json(KLINE_URL.format(symbol=symbol)).get("data") or {}
    times, closes, vols = d.get("time") or [], d.get("close") or [], d.get("vol") or []
    out = []
    for i in range(min(len(times), len(closes), len(vols))):
        close = float(closes[i])
        if close <= 0:
            continue
        out.append((int(times[i]), close, close * float(vols[i])))
    return out


def backfill(symbols: list[str] | None = None, sleep: float = 0.08) -> None:
    con = connect()
    syms = symbols if symbols is not None else universe()
    ok = failed = rows = 0
    for i, sym in enumerate(syms, 1):
        try:
            bars = fetch_bars(sym)
        except Exception:
            failed += 1
            bars = []
        if bars:
            con.executemany(
                "INSERT OR REPLACE INTO bars (symbol, ts, close, quote_vol) VALUES (?,?,?,?)",
                [(sym, t, c, q) for t, c, q in bars],
            )
            con.commit()
            ok += 1
            rows += len(bars)
        if i % 50 == 0 or i == len(syms):
            print(f"  {i}/{len(syms)}  ok={ok} failed={failed} rows={rows}", flush=True)
        time.sleep(sleep)
    con.close()
    print(f"backfill done: {ok} symbols, {rows} rows, {failed} failed")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "backfill":
        backfill()
    else:
        print(__doc__)


# ── Fuel baskets ──────────────────────────────────────────────────────────────
# Curated, not pattern-matched. Every entry was checked two ways: the lattice's
# own asset class, and the traded price against the real instrument's level.
# That second check is what matters. CVX_USDT prices at 1.49, so it is Convex
# Finance and not Chevron; Chevron is CVXSTOCK_USDT at 210. LIT_USDT prices at
# 2.30 against a lithium ETF near 50, so it is Litentry and is excluded. XLU_USDT
# and URNM_USDT did not reconcile either and are left out rather than guessed.
#
# Only 4 of the 15 WRI fuels have a proxy that survives this. Coal, Hydro, Solar,
# Biomass, Waste, Geothermal, Cogeneration, Petcoke, Storage and Wave have none
# in this universe, so they get no correlation band at all.
FUEL_TICKERS: dict[str, list[str]] = {
    "Nuclear": ["CCJSTOCK_USDT", "OKLOSTOCK_USDT", "SMRSTOCK_USDT", "VSTSTOCK_USDT"],
    "Oil": ["USOIL_USDT", "UKOIL_USDT", "XOMSTOCK_USDT", "CVXSTOCK_USDT", "OXYSTOCK_USDT"],
    "Gas": ["NGAS_USDT"],
    "Wind": ["GEVSTOCK_USDT"],
}

MIN_OVERLAP = 30


def _closes(con: sqlite3.Connection, symbol: str) -> dict[int, float]:
    return {ts: c for ts, c in con.execute(
        "SELECT ts, close FROM bars WHERE symbol=? ORDER BY ts", (symbol,))}


def basket_returns(con: sqlite3.Connection, symbols: list[str]) -> dict[int, float]:
    """Equal-weighted mean daily log return across a basket, keyed by bar time."""
    import math
    per: dict[int, list[float]] = {}
    for sym in symbols:
        closes = _closes(con, sym)
        ts = sorted(closes)
        for a, b in zip(ts, ts[1:]):
            if closes[a] > 0 and closes[b] > 0:
                per.setdefault(b, []).append(math.log(closes[b] / closes[a]))
    return {t: sum(v) / len(v) for t, v in per.items()}


def _pearson(a: list[float], b: list[float]) -> float | None:
    import math
    n = len(a)
    if n < 3:
        return None
    ma, mb = sum(a) / n, sum(b) / n
    va = sum((x - ma) ** 2 for x in a)
    vb = sum((y - mb) ** 2 for y in b)
    if va <= 0 or vb <= 0:
        return None
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / math.sqrt(va * vb)


def _pvalue(r: float, n: int) -> float:
    """Two-tailed p for a Pearson r, via the regularized incomplete beta."""
    import math
    df = n - 2
    if df <= 0 or abs(r) >= 1:
        return 0.0
    t = abs(r) * math.sqrt(df / (1 - r * r))
    x, a, b = df / (df + t * t), df / 2.0, 0.5

    def betacf(a, b, x):
        eps, tiny = 3e-16, 1e-300
        qab, qap, qam = a + b, a + 1.0, a - 1.0
        c, d = 1.0, 1.0 - qab * x / qap
        d = 1.0 / (d if abs(d) > tiny else tiny)
        h = d
        for m in range(1, 201):
            m2 = 2 * m
            for aa in (m * (b - m) * x / ((qam + m2) * (a + m2)),
                       -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))):
                d = 1.0 + aa * d
                c = 1.0 + aa / c
                d = 1.0 / (d if abs(d) > tiny else tiny)
                c = c if abs(c) > tiny else tiny
                h *= d * c
            if abs(d * c - 1.0) < eps:
                break
        return h

    lb = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    if x < (a + 1) / (a + b + 2):
        p = math.exp(a * math.log(x) + b * math.log(1 - x) - lb) * betacf(a, b, x) / a
    else:
        p = 1 - math.exp(b * math.log(1 - x) + a * math.log(x) - lb) * betacf(b, a, 1 - x) / b
    return max(0.0, min(1.0, p))


MARKET_PROXY = "QQQSTOCK_USDT"


def fuel_correlations(q: float = 0.05) -> list[dict]:
    """Pairwise fuel-basket correlations, controlled for the market, BH at q.

    Both controls are load-bearing and were added because the naive version was
    wrong in two separate ways.

    Returns, not levels: the globe's own generation series is a set of smooth
    annual trends, and correlating those measures the trends rather than any
    co-movement. Daily log returns are already differenced.

    Partial, not raw: Nuclear and Wind are equity baskets carrying market betas
    of +0.66 and +0.60, so their raw correlation of +0.53 is mostly the index.
    Controlling for QQQ drops it to +0.17. Oil<->Wind (-0.29) and Nuclear<->Oil
    (-0.28) vanish outright. Only Gas<->Oil survives, because Gas has almost no
    market beta and its link to crude is physical.
    """
    con = connect()
    series = {f: basket_returns(con, syms) for f, syms in FUEL_TICKERS.items()}
    market_ret = basket_returns(con, [MARKET_PROXY])
    con.close()
    if not market_ret:
        return []

    scored = []
    for i, fa in enumerate(sorted(series)):
        for fb in sorted(series)[i + 1:]:
            shared = sorted(set(series[fa]) & set(series[fb]) & set(market_ret))
            n = len(shared)
            if n < MIN_OVERLAP:
                continue
            a = [series[fa][t] for t in shared]
            b = [series[fb][t] for t in shared]
            z = [market_ret[t] for t in shared]
            rab, raz, rbz = _pearson(a, b), _pearson(a, z), _pearson(b, z)
            if None in (rab, raz, rbz):
                continue
            denom = ((1 - raz ** 2) * (1 - rbz ** 2)) ** 0.5
            if denom <= 0:
                continue
            partial = (rab - raz * rbz) / denom
            scored.append({
                "a": fa, "b": fb,
                "r": round(partial, 4), "r_raw": round(rab, 4),
                "n": n, "p": _pvalue(partial, n - 1),  # one control costs a df
            })

    scored.sort(key=lambda d: d["p"])
    m = len(scored)
    keep = 0
    for i, d in enumerate(scored, 1):
        if d["p"] <= q * i / m:
            keep = i
    for i, d in enumerate(scored, 1):
        d["significant"] = i <= keep
        d["p"] = float(f"{d['p']:.3g}")
    return scored
