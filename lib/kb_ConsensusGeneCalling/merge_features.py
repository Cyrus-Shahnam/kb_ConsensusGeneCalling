# -*- coding: utf-8 -*-
from __future__ import annotations

from typing import Dict, List, Tuple

TYPE_PRIORITY = {
    "CDS": 0,
    "tRNA": 1,
    "rRNA": 2,
    "ncRNA": 3,
}

def _norm_interval(f: Dict) -> Tuple[str, int, int]:
    """
    Always treat start<=stop in comparisons; caller outputs can be mixed.
    """
    s = int(f["start"])
    e = int(f["stop"])
    if s <= e:
        return f["contig"], s, e
    return f["contig"], e, s

def _contig_rank(contig_order: List[str]) -> Dict[str, int]:
    return {c: i for i, c in enumerate(contig_order)}

def merge_overlapping_features(
    features: List[Dict],
    contig_order: List[str],
    coordinate_system: str = "1-based_closed",
) -> List[Dict]:
    """
    Collapse overlapping features into "regions".
    - If CDS/tRNA/rRNA overlap at same locus -> region.type becomes a joined label (e.g., "CDS+tRNA+rRNA")
    - Ordering: contig order, then start, then priority (CDS first), then stop.
    - Keeps caller provenance via a "callers" field and "types" field.
    """

    if not features:
        return []

    rank = _contig_rank(contig_order)

    # sort features by locus, then by your type priority (CDS>tRNA>rRNA)
    def _sort_key(f):
        contig, s, e = _norm_interval(f)
        tp = TYPE_PRIORITY.get(f.get("type", "ncRNA"), 99)
        return (rank.get(contig, 10**9), contig, s, tp, e)

    feats = sorted(features, key=_sort_key)

    merged: List[Dict] = []
    cur = None  # current region dict

    def _start_new_region(f: Dict) -> Dict:
        contig, s, e = _norm_interval(f)
        region = {
            "contig": contig,
            "start": s,
            "stop": e,
            # keep strands if you want; for mixed strands just mark "."
            "direction": f.get("direction", "."),
            "types": set([f.get("type", "ncRNA")]),
            "callers": set([f.get("caller", "")]) if f.get("caller") else set(),
            # store originals if you ever want debugging
            "_members": [f],
        }
        return region

    def _can_merge(region: Dict, f: Dict) -> bool:
        contig, s, e = _norm_interval(f)
        if contig != region["contig"]:
            return False
        # overlap if s <= region.stop AND e >= region.start
        return s <= int(region["stop"]) and e >= int(region["start"])

    for f in feats:
        if cur is None:
            cur = _start_new_region(f)
            continue

        if _can_merge(cur, f):
            contig, s, e = _norm_interval(f)
            cur["start"] = min(int(cur["start"]), s)
            cur["stop"] = max(int(cur["stop"]), e)

            # direction handling
            d0 = cur.get("direction", ".")
            d1 = f.get("direction", ".")
            if d0 == ".":
                cur["direction"] = d1
            elif d1 != "." and d1 != d0:
                cur["direction"] = "."

            cur["types"].add(f.get("type", "ncRNA"))
            if f.get("caller"):
                cur["callers"].add(f["caller"])
            cur["_members"].append(f)
        else:
            merged.append(cur)
            cur = _start_new_region(f)

    if cur is not None:
        merged.append(cur)

    # finalize region objects to match your output schema
    out: List[Dict] = []
    for i, r in enumerate(merged, start=1):
        types_sorted = sorted(list(r["types"]), key=lambda t: TYPE_PRIORITY.get(t, 99))
        callers_sorted = sorted([c for c in r["callers"] if c])

        out.append({
            "caller": ",".join(callers_sorted) if callers_sorted else "",
            "ID": f"region.{i}",
            "contig": r["contig"],
            "start": int(r["start"]),
            "stop": int(r["stop"]),
            "direction": r.get("direction", "."),
            "DNA sequence": "",  # regions don't carry sequence
            "type": "+".join(types_sorted),
        })

    return out
