"""Component-overlap histories for the supra-threshold region.

A per-frame component COUNT cannot distinguish ten stable regions from ten regions
churning every frame. "Consolidates to 10 and holds for 85 frames" is a claim about
a count; whether it is also a claim about identity is exactly what this measures.

For each consecutive pair of frames, intersect every component's voxel set with every
component's voxel set in the next frame. A component inherits the lineage of whichever
predecessor it shares the most voxels with; ties and contested inheritances are resolved
largest-claim-first, and everything unclaimed is a birth.

Both adjacency rules are reported, always labelled. ParaView's Connectivity filter joins
voxels touching at a corner (26-connectivity); scipy's default structure is 6. The same
region gives different counts under each, so a bare component count is not a number.

usage: component_overlap.py <inert_signal_dir> [--out <dir>]
"""
import sys, os, json
import numpy as np
from scipy import ndimage
from vti_read import read_vti

S6 = ndimage.generate_binary_structure(3, 1)   # face-sharing, 6-connectivity
S26 = ndimage.generate_binary_structure(3, 3)  # point-sharing, 26-connectivity, = ParaView


def label_frames(root, struct, n=101):
    """Yield (mcs, labelled array, component count) per frame."""
    for t in range(n):
        a, _, _ = read_vti(os.path.join(root, "paraview", "signal_mcs%06d.vti" % t))
        lab, k = ndimage.label(a["occupied_above_threshold"] != 0, structure=struct)
        yield t, lab, k


def overlap_history(root, struct, n=101):
    counts, events, lineage_of_prev, next_lineage = [], [], {}, 0
    lineage_life = {}     # lineage id -> [first_mcs, last_mcs, frames_seen]
    prev_lab = None

    for t, lab, k in label_frames(root, struct, n):
        counts.append(k)
        cur_sizes = {i: int((lab == i).sum()) for i in range(1, k + 1)}

        if prev_lab is None or k == 0:
            lineage_of_cur = {}
            for i in range(1, k + 1):
                lineage_of_cur[i] = next_lineage
                lineage_life[next_lineage] = [t, t, 1]
                next_lineage += 1
            events.append({"mcs": t, "births": k, "deaths": 0, "merges": 0, "splits": 0,
                           "carried": 0})
        else:
            # intersection counts, only where both are non-zero
            both = (prev_lab > 0) & (lab > 0)
            pairs = {}
            if both.any():
                pl, cl = prev_lab[both], lab[both]
                for p, c in zip(pl.ravel(), cl.ravel()):
                    pairs[(int(p), int(c))] = pairs.get((int(p), int(c)), 0) + 1

            # largest claim wins; a predecessor's lineage is inherited once
            claims = sorted(pairs.items(), key=lambda kv: -kv[1])
            lineage_of_cur, used_pred, taken_cur = {}, set(), set()
            for (p, c), _ in claims:
                if p in used_pred or c in taken_cur:
                    continue
                lineage_of_cur[c] = lineage_of_prev[p]
                used_pred.add(p); taken_cur.add(c)

            births = 0
            for i in range(1, k + 1):
                if i not in lineage_of_cur:
                    lineage_of_cur[i] = next_lineage
                    lineage_life[next_lineage] = [t, t, 0]
                    next_lineage += 1
                    births += 1

            preds_with_succ = {p for (p, c) in pairs}
            succs_per_pred = {}
            preds_per_succ = {}
            for (p, c) in pairs:
                succs_per_pred.setdefault(p, set()).add(c)
                preds_per_succ.setdefault(c, set()).add(p)
            deaths = len(set(lineage_of_prev)) - len(preds_with_succ)
            merges = sum(1 for c, ps in preds_per_succ.items() if len(ps) > 1)
            splits = sum(1 for p, cs in succs_per_pred.items() if len(cs) > 1)
            events.append({"mcs": t, "births": births, "deaths": deaths,
                           "merges": merges, "splits": splits,
                           "carried": len(used_pred)})

        for i, lin in lineage_of_cur.items():
            lineage_life[lin][1] = t
            lineage_life[lin][2] += 1
        lineage_of_prev = lineage_of_cur
        prev_lab = lab

    return counts, events, lineage_life


def main(argv):
    root = os.path.abspath(argv[1])
    out = argv[argv.index("--out") + 1] if "--out" in argv else os.path.dirname(os.path.abspath(__file__))
    report = {"evidence_root": root, "adjacency": {}}

    for name, struct in (("26-connectivity (point-sharing; ParaView Connectivity)", S26),
                         ("6-connectivity (face-sharing; scipy default)", S6)):
        counts, events, life = overlap_history(root, struct)
        tail = [e for e in events if e["mcs"] >= 15]
        churn = sum(e["births"] + e["deaths"] for e in tail)
        alive_at_100 = [l for l, (f, la, _) in life.items() if la == 100]
        born_before_15 = [l for l in alive_at_100 if life[l][0] < 15]
        report["adjacency"][name] = {
            "counts_by_mcs": counts,
            "count_at": {str(m): counts[m] for m in (7, 8, 9, 11, 15, 30, 50, 100)},
            "peak": {"mcs": int(np.argmax(counts)), "count": int(max(counts))},
            "flat_from_15": len(set(counts[15:])) == 1,
            "distinct_lineages_ever": len(life),
            "lineages_alive_at_mcs100": len(alive_at_100),
            "of_those_born_before_mcs15": len(born_before_15),
            "births_plus_deaths_after_mcs15": churn,
            "events": events,
            "lineage_lifespans": {str(l): {"first": f, "last": la, "frames": n}
                                  for l, (f, la, n) in sorted(life.items())},
        }
        print("\n=== %s ===" % name)
        print("  counts at 7/8/9/11/15/30/50/100 : %s" %
              [counts[m] for m in (7, 8, 9, 11, 15, 30, 50, 100)])
        print("  peak                            : %d at MCS %d" % (max(counts), int(np.argmax(counts))))
        print("  count flat from MCS 15 on       : %s" % (len(set(counts[15:])) == 1))
        print("  distinct lineages ever          : %d" % len(life))
        print("  alive at MCS 100                : %d" % len(alive_at_100))
        print("    of those, born before MCS 15  : %d" % len(born_before_15))
        print("  births+deaths after MCS 15      : %d   <-- churn behind a flat count" % churn)

    dest = os.path.join(out, "component_overlap.json")
    json.dump(report, open(dest, "w"), indent=2)
    print("\nwrote %s" % dest)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
