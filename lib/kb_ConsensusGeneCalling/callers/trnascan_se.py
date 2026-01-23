# -*- coding: utf-8 -*-
from __future__ import annotations

import os
from typing import List, Dict, Optional

from kb_ConsensusGeneCalling.gene_calling_io import ensure_dir, run_cmd


def run_trnascan_se(
    fasta_path: str,
    out_dir: str,
    mode: str = "bacterial",
) -> Dict[str, str]:
    """
    Run tRNAscan-SE and write tabular output.
    Output:
      - trnascan.tsv
    """
    ensure_dir(out_dir)
    out_tsv = os.path.join(out_dir, "trnascan.tsv")

    # tRNAscan-SE 2.x common flags
    mode_flag = "-B"
    m = (mode or "").lower().strip()
    if m.startswith("arch"):
        mode_flag = "-A"
    elif m.startswith("euk"):
        mode_flag = "-E"

    cmd = [
        "/opt/conda3/bin/conda", "run", "-n", "gene_calling",
        "tRNAscan-SE",
        mode_flag,
        "-o", out_tsv,
        fasta_path,
    ]
    run_cmd(cmd)
    return {"tsv": out_tsv}


def _try_int(x: str) -> Optional[int]:
    try:
        return int(x)
    except Exception:
        return None


def parse_trnascan_tsv(tsv_path: str) -> List[Dict]:
    """
    Parses TSV like:

      # tRNAscan-SE test example
      contigA  1  10  50  tRNA-Leu  TAA  TA
      contigB  2  120 80  tRNA-Gly  CCC  CC

    Returns raw feats with:
      contig, start, stop, strand, type="tRNA", id="contigX.tRNA.N"
    Coordinates are 1-based inclusive.
    """
    feats: List[Dict] = []
    if not os.path.exists(tsv_path):
        return feats

    with open(tsv_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # tolerate tabs OR spaces
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) < 4:
                continue

            contig = parts[0]
            hit_num = parts[1] if len(parts) >= 2 else ""
            begin = _try_int(parts[2]) if len(parts) >= 3 else None
            end = _try_int(parts[3]) if len(parts) >= 4 else None
            if begin is None or end is None:
                continue

            if begin <= end:
                start, stop, strand = begin, end, "+"
            else:
                start, stop, strand = end, begin, "-"

            feats.append({
                "contig": contig,
                "start": int(start),
                "stop": int(stop),
                "strand": strand,
                "type": "tRNA",
                "_hit_num": str(hit_num),
                "source": "tRNAscan-SE",
            })

    # deterministic ordering
    feats.sort(key=lambda r: (r["contig"], int(r["start"]), int(r["stop"])))

    # assign IDs exactly as tests expect: contigA.tRNA.1, contigB.tRNA.1, ...
    per_contig_counter: Dict[str, int] = {}
    for rf in feats:
        c = rf["contig"]
        per_contig_counter[c] = per_contig_counter.get(c, 0) + 1
        rf["id"] = f"{c}.tRNA.{per_contig_counter[c]}"
        rf.pop("_hit_num", None)

    return feats


# Compatibility alias (if anything imports this older name)
def parse_trnascan_output(path: str) -> List[Dict]:
    return parse_trnascan_tsv(path)
