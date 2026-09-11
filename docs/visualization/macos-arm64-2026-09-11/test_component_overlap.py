#!/usr/bin/env python3
"""Fixtures for the overlap accounting. Synthetic frames, no VTI, no evidence bundle.

Every one of these failed against the pre-repair implementation. They exist because the
production numbers looked healthy while three separate counters were wrong: a lineage
present for one frame reported two, an all-vanish frame reported zero deaths, and every
lineage that retired by losing a merge claim was counted nowhere at all.
"""
import sys
import numpy as np
from component_overlap import overlap_history

FAILS = []


def check(name, cond, got):
    print(("  PASS  " if cond else "  FAIL  ") + name + "  ->  " + str(got))
    if not cond:
        FAILS.append(name)


def frame(mcs, blocks, shape=(9, 1, 1)):
    """blocks: list of (start, stop) index ranges, each becoming one labelled component."""
    lab = np.zeros(shape, dtype=int)
    for i, (a, b) in enumerate(blocks, start=1):
        lab[a:b, 0, 0] = i
    return (mcs, lab, len(blocks))


def run(frames):
    return overlap_history(None, None, frames=frames)


print("=== one component, one frame: lifespan is 1, not 2 ===")
counts, events, life = run([frame(0, [(0, 2)])])
check("frames_seen == 1", life[0][2] == 1, life[0][2])
check("counts == [1]", counts == [1], counts)

print("\n=== one component, then an empty frame: one disappearance ===")
counts, events, life = run([frame(0, [(0, 2)]), frame(1, [])])
e = events[1]
check("deaths == 1", e["deaths"] == 1, e["deaths"])
check("retired_by_merge == 0", e["retired_by_merge"] == 0, e["retired_by_merge"])
check("carried == 0", e["carried"] == 0, e["carried"])

print("\n=== two components joining into one: 1 merge, 1 retirement, 0 deaths ===")
counts, events, life = run([frame(0, [(0, 2), (4, 6)]), frame(1, [(0, 6)])])
e = events[1]
check("merges == 1", e["merges"] == 1, e["merges"])
check("retired_by_merge == 1", e["retired_by_merge"] == 1, e["retired_by_merge"])
check("deaths == 0 (neither disappeared)", e["deaths"] == 0, e["deaths"])
check("carried == 1", e["carried"] == 1, e["carried"])
check("births == 0", e["births"] == 0, e["births"])

print("\n=== one component splitting into two: 1 split, 1 birth, 0 deaths ===")
counts, events, life = run([frame(0, [(0, 6)]), frame(1, [(0, 2), (4, 6)])])
e = events[1]
check("splits == 1", e["splits"] == 1, e["splits"])
check("births == 1", e["births"] == 1, e["births"])
check("deaths == 0", e["deaths"] == 0, e["deaths"])
check("retired_by_merge == 0", e["retired_by_merge"] == 0, e["retired_by_merge"])

print("\n=== competing overlap: larger claim wins, loser retires ===")
# prev: A=[0,5) B=[6,8);  cur: one component [0,8) overlapping A by 5 and B by 2
counts, events, life = run([frame(0, [(0, 5), (6, 8)]), frame(1, [(0, 8)])])
e = events[1]
check("carried == 1", e["carried"] == 1, e["carried"])
check("retired_by_merge == 1", e["retired_by_merge"] == 1, e["retired_by_merge"])
check("surviving lineage is the LARGER predecessor (id 0)",
      life[0][1] == 1 and life[1][1] == 0, "lineage0.last=%d lineage1.last=%d" % (life[0][1], life[1][1]))

print("\n=== three into one: 1 merge, 2 retirements ===")
counts, events, life = run([frame(0, [(0, 2), (3, 5), (6, 8)]), frame(1, [(0, 8)])])
e = events[1]
check("merges == 1", e["merges"] == 1, e["merges"])
check("retired_by_merge == 2", e["retired_by_merge"] == 2, e["retired_by_merge"])

print("\n=== disappearance and merge in the same transition ===")
# prev: A=[0,2) B=[3,5) C=[7,9);  cur: [0,5) merges A+B, C vanishes
counts, events, life = run([frame(0, [(0, 2), (3, 5), (7, 9)]), frame(1, [(0, 5)])])
e = events[1]
check("deaths == 1 (C disappeared)", e["deaths"] == 1, e["deaths"])
check("retired_by_merge == 1 (B lost)", e["retired_by_merge"] == 1, e["retired_by_merge"])
check("carried == 1", e["carried"] == 1, e["carried"])

print("\n%d checks, %d failures" % (7 + 2 + 5 + 4 + 3 + 2 + 3 - len(FAILS) + len(FAILS), len(FAILS)))
sys.exit(1 if FAILS else 0)
