# -*- coding: utf-8 -*-
from __future__ import annotations

import os
from typing import List, Dict

from kb_ConsensusGeneCalling.gene_calling_io import ensure_dir, run_cmd


def run_trnascan_se(fasta_path: str, out_dir: str) -> dict:
    """
    Runs tRNAscan-SE and writes outputs into out_dir.

    NOTE: CLI flags vary slightly across versions. This is a reasonable default for v2.x,
    but if your environment differs, adjust flags accordingly.
    """
    ensure_dir(out_dir)
    gff_path = os.path.join(out_dir, "trnascan.gff")
    txt_path = os.path.join(out_dir, "trnascan.txt")

    # Prefer producing a GFF file if supported.
    # Some builds use '--output' and '--gff'; others use '-o' and '--gff'.
    cmd = [
        "/opt/conda3/bin/conda", "run", "-n", "gene_calling",
        "tRNAscan-SE",
        "--gff", gff_path,
        "-o", txt_path,
        fasta_path
    ]
    run_cmd(cmd)

    return {"gff": gff_path, "txt": txt_path}


def parse_trnascan_gff(gff_path: str) -> List[Dict]:
    """
    Parse GFF-like output into raw feats suitable for normalization later.
    """
    feats: List[Dict] = []
    with open(gff_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                continue
            contig, source, ftype, start, stop, score, strand, phase, attrs = parts

            # Standardize type
            out_type = "tRNA" if "trna" in ftype.lower() else "ncRNA"

            attr_map = {}
            for tok in attrs.split(";"):
                tok = tok.strip()
                if not tok:
                    continue
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    attr_map[k] = v

            feats.append({
                "contig": contig,
                "start": int(start),
                "stop": int(stop),
                "strand": strand,
                "type": out_type,
                "id": attr_map.get("ID") or attr_map.get("Name"),
                "attributes": attr_map
            })
    return feats
