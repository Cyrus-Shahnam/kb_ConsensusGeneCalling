import os
from typing import Dict, List, Tuple

from installed_clients.ServiceWizardClient import ServiceWizard
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil


def _get_service_url(callback_url: str, module_name: str) -> str:
    sw = ServiceWizard(callback_url)
    return sw.get_service_status({"module_name": module_name})["url"]


def _read_fasta_records(fasta_path: str):
    """
    Minimal FASTA reader (no external deps).
    Yields (header, seq) where header is up to first whitespace.
    """
    header = None
    seq_chunks = []
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def _load_contigs_from_fasta(
    fasta_path: str,
    min_contig_length: int = 0,
) -> Tuple[List[str], Dict[str, str]]:
    """
    Returns:
      contig_order: list of contig IDs in file order
      contigs: dict contig_id -> sequence (uppercased)
    Filters contigs shorter than min_contig_length.
    """
    contig_order: List[str] = []
    contigs: Dict[str, str] = {}

    for cid, seq in _read_fasta_records(fasta_path):
        seq = seq.upper()
        if min_contig_length and len(seq) < int(min_contig_length):
            continue
        contig_order.append(cid)
        contigs[cid] = seq

    return contig_order, contigs


def fetch_assembly_as_fasta(
    ctx,
    assembly_ref: str,
    scratch: str,
    min_contig_length: int = 0,
):
    """
    Fetch contigs for a KBase Assembly ref.

    Appdev/prod-friendly approach:
      1) Try AssemblyUtil.get_assembly_as_fasta (recommended)
      2) Fallback to DFU.get_objects for legacy/edge cases

    Returns:
      contig_order (list[str]), contigs (dict[str,str])
    """
    callback_url = os.environ.get("SDK_CALLBACK_URL")
    if not callback_url:
        raise ValueError("SDK_CALLBACK_URL is not set; cannot resolve service URLs.")

    # 1) Preferred: AssemblyUtil
    try:
        au_url = _get_service_url(callback_url, "AssemblyUtil")
        au = AssemblyUtil(au_url)
        res = au.get_assembly_as_fasta({"ref": assembly_ref})
        fasta_path = res["path"]  # local file path inside container
        contig_order, contigs = _load_contigs_from_fasta(
            fasta_path, min_contig_length=min_contig_length
        )
        if not contig_order:
            raise ValueError(
                f"AssemblyUtil returned FASTA but all contigs were filtered out "
                f"(min_contig_length={min_contig_length})."
            )
        return contig_order, contigs
    except Exception as e:
        # continue to DFU fallback
        last_err = e

    # 2) Fallback: DFU get_objects (less reliable for assemblies)
    try:
        dfu = DataFileUtil(callback_url)
        obj = dfu.get_objects({"object_refs": [assembly_ref]})["data"][0]["data"]

        # Assembly objects vary; try common places sequences might exist.
        # If your old DFU fallback logic was more complete, keep/merge it here.
        contigs = {}
        contig_order = []

        # Some assemblies may contain 'contigs' array with 'id' and 'sequence'
        if isinstance(obj, dict) and "contigs" in obj and isinstance(obj["contigs"], list):
            for c in obj["contigs"]:
                cid = c.get("id") or c.get("contig_id")
                seq = c.get("sequence")
                if not cid or not seq:
                    continue
                seq = str(seq).upper()
                if min_contig_length and len(seq) < int(min_contig_length):
                    continue
                contig_order.append(cid)
                contigs[cid] = seq

        if not contig_order:
            raise ValueError(
                "Could not obtain contig sequences from Assembly via DFU fallback, "
                "or all contigs were filtered out. If using KBase Assembly objects, "
                "AssemblyUtil is recommended."
            )

        return contig_order, contigs

    except Exception:
        # show the AssemblyUtil error too — it’s usually the real clue
        raise ValueError(
            "Could not obtain contig sequences from Assembly via AssemblyUtil and DFU fallback. "
            f"AssemblyUtil error was: {repr(last_err)}"
        )
