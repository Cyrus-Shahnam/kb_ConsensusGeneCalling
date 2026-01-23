# -*- coding: utf-8 -*-
from __future__ import annotations

import json
from typing import List, Dict

FEATURE_TSV_COLUMNS = [
    "caller",          # top caller by your priority
    "callers",         # all callers that overlapped in this region
    "ID",
    "contig",
    "start",
    "stop",
    "direction",
    "type",            # e.g. CDS+tRNA+rRNA
    "DNA sequence",
]

def _clean(v) -> str:
    if v is None:
        return ""
    s = str(v)
    return s.replace("\t", " ").replace("\n", "").replace("\r", "")

def write_features_tsv(features: List[Dict], out_tsv: str) -> None:
    with open(out_tsv, "w") as f:
        f.write("\t".join(FEATURE_TSV_COLUMNS) + "\n")
        for feat in features:
            row = [_clean(feat.get(col, "")) for col in FEATURE_TSV_COLUMNS]
            f.write("\t".join(row) + "\n")

def write_features_json(features: List[Dict], out_json: str) -> None:
    with open(out_json, "w") as f:
        json.dump(features, f, indent=2, sort_keys=False)
