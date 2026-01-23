# -*- coding: utf-8 -*-
from __future__ import annotations

import os
import shutil
import subprocess
from typing import List, Dict, Tuple, Optional


def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)


def write_fasta(contigs: List[Dict[str, str]], out_fa: str) -> List[str]:
    """
    contigs: list of dicts: [{"id": str, "seq": str}, ...]
    Returns contig id order.
    """
    ids: List[str] = []
    with open(out_fa, "w") as f:
        for c in contigs:
            cid = c["id"]
            seq = c["seq"].replace(" ", "").replace("\n", "").upper()
            ids.append(cid)
            f.write(f">{cid}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")
    return ids


def parse_fasta_ids(in_fa: str) -> List[str]:
    ids: List[str] = []
    with open(in_fa) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def stage_input_fasta(params: dict, scratch: str) -> Tuple[str, List[str]]:
    """
    MVP path: if params["input_fasta_path"] is provided, copy it into scratch and return (fasta_path, contig_order).
    """
    src = params.get("input_fasta_path")
    if not src:
        raise ValueError("Missing input_fasta_path for local mode")
    if not os.path.exists(src):
        raise ValueError(f"input_fasta_path does not exist: {src}")

    ensure_dir(scratch)
    dst = os.path.join(scratch, "input.fna")
    shutil.copyfile(src, dst)
    order = parse_fasta_ids(dst)
    if not order:
        raise ValueError("No contigs found in input FASTA")
    return dst, order


def _filter_fasta_by_min_len(in_fa: str, out_fa: str, min_len: int) -> List[str]:
    """
    Stream FASTA, keep only sequences with length >= min_len.
    Returns kept contig order.
    """
    kept_order: List[str] = []

    def flush(cur_id: Optional[str], buf: List[str], out_handle) -> None:
        if cur_id is None:
            return
        seq = "".join(buf).replace(" ", "").replace("\n", "").upper()
        if len(seq) >= min_len:
            kept_order.append(cur_id)
            out_handle.write(f">{cur_id}\n")
            for i in range(0, len(seq), 60):
                out_handle.write(seq[i:i + 60] + "\n")

    with open(in_fa) as fin, open(out_fa, "w") as fout:
        cur_id: Optional[str] = None
        buf: List[str] = []
        for line in fin:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush(cur_id, buf, fout)
                cur_id = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        flush(cur_id, buf, fout)

    return kept_order


def fetch_assembly_as_fasta(
    dfu,
    assembly_ref: str,
    out_fa: Optional[str] = None,
    out_fasta_path: Optional[str] = None,
    min_contig_length: int = 0
) -> List[str]:
    """
    Fetch assembly FASTA and optionally filter contigs by min_contig_length.

    Supports BOTH call styles:
      - fetch_assembly_as_fasta(dfu=..., assembly_ref=..., out_fa="x.fna")
      - fetch_assembly_as_fasta(dfu=..., assembly_ref=..., out_fasta_path="x.fna", min_contig_length=200)

    Returns contig id order (after filtering).
    """
    # accept either param name
    out_path = out_fasta_path or out_fa
    if not out_path:
        raise ValueError("Missing out_fa/out_fasta_path")

    # 1) Try AssemblyUtil if available
    tmp_path = out_path
    try:
        from installed_clients.AssemblyUtilClient import AssemblyUtil  # type: ignore

        au = AssemblyUtil(os.environ["SDK_CALLBACK_URL"])
        # AssemblyUtil writes fasta to provided filename
        au.get_assembly_as_fasta({"ref": assembly_ref, "filename": tmp_path})

        if min_contig_length and min_contig_length > 0:
            # Filter into a temp file then replace
            filtered = tmp_path + ".filtered"
            kept = _filter_fasta_by_min_len(tmp_path, filtered, min_contig_length)
            shutil.move(filtered, tmp_path)
            if not kept:
                raise ValueError(
                    f"No contigs >= min_contig_length={min_contig_length} "
                    f"after filtering AssemblyUtil FASTA"
                )
            return kept

        order = parse_fasta_ids(tmp_path)
        if not order:
            raise ValueError("No contigs found in AssemblyUtil FASTA")
        return order

    except Exception:
        # fall through to DFU inline fallback
        pass

    # 2) DFU fallback (only works if contigs include inline sequences)
    obj = dfu.get_objects({"object_refs": [assembly_ref]})["data"][0]
    data = obj["data"]

    contigs: List[Dict[str, str]] = []
    if "contigs" in data and isinstance(data["contigs"], list) and data["contigs"]:
        for c in data["contigs"]:
            seq = c.get("sequence")
            if not seq:
                continue
            if min_contig_length and len(seq) < int(min_contig_length):
                continue

            cid = c.get("id") or c.get("contig_id") or c.get("name")
            if not cid:
                cid = f"contig_{len(contigs) + 1}"
            contigs.append({"id": cid, "seq": seq})

    if not contigs:
        raise ValueError(
            "Could not obtain contig sequences from Assembly via DFU fallback, "
            "or all contigs were filtered out. "
            "If using KBase Assembly objects, AssemblyUtil is recommended."
        )

    return write_fasta(contigs, tmp_path)


def run_cmd(cmd: List[str], cwd: Optional[str] = None) -> None:
    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            f"  cmd: {' '.join(cmd)}\n"
            f"  exit: {p.returncode}\n"
            f"  stdout:\n{p.stdout}\n"
            f"  stderr:\n{p.stderr}\n"
        )
