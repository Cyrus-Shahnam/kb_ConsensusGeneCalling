# -*- coding: utf-8 -*-
from __future__ import annotations

import os
from typing import List, Dict

from kb_ConsensusGeneCalling.gene_calling_io import ensure_dir, run_cmd


def run_barrnap(
    fasta_path: str,
    out_dir: str,
    kingdom: str = "bac",
) -> Dict[str, str]:
    """
    Run barrnap and write output GFF.
    Output:
      - barrnap.gff
    """
    ensure_dir(out_dir)
    out_gff = os.path.join(out_dir, "barrnap.gff")

    k = (kingdom or "").lower().strip()
    if k not in ("bac", "arc", "euk", "mito"):
        k = "bac"

    # barrnap emits GFF to stdout
    cmd = f"/opt/conda3/bin/conda run -n gene_calling barrnap --kingdom {k} {fasta_path} > {out_gff}"
    run_cmd(["/bin/bash", "-lc", cmd])

    return {"gff": out_gff}


def parse_barrnap_gff(gff_path: str) -> List[Dict]:
    """
    Parses barrnap GFF like:

      ##gff-version 3
      contigA  barrnap  16S_rRNA  5   120  .  +  .  ID=rrna1;Name=16S_rRNA
      contigB  barrnap  rRNA      10  90   .  -  .  ID=rrna2;Name=23S_rRNA

    Returns raw feats with:
      contig, start, stop, strand, type="rRNA", id
    Coordinates are 1-based inclusive.
    """
    feats: List[Dict] = []
    if not os.path.exists(gff_path):
        return feats

    with open(gff_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                parts = line.split()
            if len(parts) < 9:
                continue

            contig, source, ftype, start, stop, score, strand, phase, attrs = parts[:9]

            # barrnap feature types include: 16S_rRNA, 23S_rRNA, 5S_rRNA, rRNA
            if "rrna" not in (ftype or "").lower():
                continue

            try:
                start_i = int(start)
                stop_i = int(stop)
            except Exception:
                continue

            attr_map: Dict[str, str] = {}
            for tok in (attrs or "").split(";"):
                tok = tok.strip()
                if not tok:
                    continue
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    attr_map[k] = v

            feat_id = attr_map.get("ID") or f"rrna_{contig}_{start_i}_{stop_i}"

            feats.append({
                "contig": contig,
                "start": start_i,
                "stop": stop_i,
                "strand": "+" if strand == "+" else "-",
                "type": "rRNA",
                "id": feat_id,
                "source": source or "barrnap",
            })

    feats.sort(key=lambda r: (r["contig"], int(r["start"]), int(r["stop"])))
    return feats
