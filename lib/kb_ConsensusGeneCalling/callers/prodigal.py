# -*- coding: utf-8 -*-
from __future__ import annotations

import os
from typing import Dict, List, Optional

from kb_ConsensusGeneCalling.gene_calling_io import run_cmd, ensure_dir


def _fasta_total_and_max_len(fasta_path: str) -> tuple[int, int]:
    total = 0
    max_len = 0
    cur = 0
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur > 0:
                    total += cur
                    max_len = max(max_len, cur)
                cur = 0
            else:
                cur += len(line)
        if cur > 0:
            total += cur
            max_len = max(max_len, cur)
    return total, max_len


def run_prodigal(
    fasta_path: str,
    out_dir: str,
    genetic_code: int = 11,
    mode: str = "auto"  # "single" or "meta" or "auto"
) -> dict:
    """
    Produces:
      - genes.gff (GFF)
      - proteins.faa (FAA)
      - genes.ffn (FFN nucleotide sequences on coding strand)
    """
    ensure_dir(out_dir)

    gff = os.path.join(out_dir, "genes.gff")
    faa = os.path.join(out_dir, "proteins.faa")
    ffn = os.path.join(out_dir, "genes.ffn")

    if mode not in ("single", "meta", "auto"):
        raise ValueError(f"Invalid prodigal mode: {mode}")

    chosen_mode = mode
    if mode == "auto":
        total_len, max_len = _fasta_total_and_max_len(fasta_path)
        # Prodigal single-genome training needs >= ~20kb.
        chosen_mode = "meta" if (total_len < 20000 or max_len < 20000) else "single"

    cmd = [
        "/opt/conda3/bin/conda", "run", "-n", "gene_calling",
        "prodigal",
        "-i", fasta_path,
        "-o", gff,
        "-a", faa,
        "-d", ffn,
        "-f", "gff",
        "-g", str(genetic_code),
        "-p", "meta" if chosen_mode == "meta" else "single",
    ]

    run_cmd(cmd)

    return {"gff": gff, "faa": faa, "ffn": ffn, "mode": chosen_mode}


def parse_fasta_map(fa_path: str) -> dict[str, str]:
    """
    Parses FASTA into {id: sequence}. ID = first token on header line.
    """
    seqs: dict[str, list[str]] = {}
    cur = None
    with open(fa_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            else:
                if cur is None:
                    raise ValueError(f"FASTA parse error in {fa_path}")
                seqs[cur].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def parse_prodigal_gff(gff_path: str) -> list[dict]:
    """
    Return raw features parsed from Prodigal GFF.
    Prodigal GFF coordinates are 1-based inclusive.
    """
    feats: list[dict] = []

    with open(gff_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) != 9:
                continue

            contig, source, ftype, start, stop, score, strand, phase, attrs = parts

            if ftype != "CDS":
                continue

            start_i = int(start)
            stop_i = int(stop)

            attr_map = {}
            for token in attrs.split(";"):
                token = token.strip()
                if not token:
                    continue
                if "=" in token:
                    k, v = token.split("=", 1)
                    attr_map[k] = v

            prod_id = attr_map.get("ID")
            partial = attr_map.get("partial")
            start_type = attr_map.get("start_type")

            if not prod_id:
                prod_id = f"{contig}_{start_i}_{stop_i}_{strand}"

            feats.append({
                "contig": contig,
                "start": start_i,
                "stop": stop_i,
                "strand": strand,
                "prodigal_id": prod_id,
                "partial": partial,
                "type": "CDS",
                "source": source,
                "start_type": start_type,
            })

    return feats
