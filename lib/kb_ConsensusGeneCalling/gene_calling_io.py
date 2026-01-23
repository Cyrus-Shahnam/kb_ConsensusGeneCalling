import os
import subprocess
from typing import Dict, List, Tuple

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil

def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path

def run_cmd(cmd, cwd=None, env=None, check=True):
    """
    Run a command (list or string). Returns stdout as text.
    Raises CalledProcessError if check=True and command fails.
    """
    shell = isinstance(cmd, str)
    p = subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        shell=shell,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if check and p.returncode != 0:
        raise subprocess.CalledProcessError(
            p.returncode, cmd, output=p.stdout, stderr=p.stderr
        )
    return p.stdout

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
    ctx=None,
    assembly_ref: str = None,
    scratch: str = None,
    min_contig_length: int = 0,
    dfu=None,
    out_fasta_path: str = None,
    **kwargs,
):

    """
    Fetch contigs for a KBase Assembly ref.

    Accepts extra kwargs for compatibility with older/newer callers:
      - dfu: optional DataFileUtil client
      - out_fasta_path: optional path to write the fetched FASTA
      - **kwargs: ignored

    NO ServiceWizard version:
      1) Try AssemblyUtil.get_assembly_as_fasta using SDK_CALLBACK_URL
      2) Fallback to DFU.get_objects (best-effort)

    Returns:
      contig_order (list[str]), contigs (dict[str,str])
    """
    callback_url = os.environ.get("SDK_CALLBACK_URL")
    if not callback_url:
        raise ValueError("SDK_CALLBACK_URL is not set; cannot call AssemblyUtil/DFU.")

    last_err = None

    # 1) Preferred: AssemblyUtil via callback URL
    try:
        au = AssemblyUtil(callback_url)
        res = au.get_assembly_as_fasta({"ref": assembly_ref})
        fasta_path = res["path"]  # local file path inside container

        # Optional: copy FASTA to requested location
        if out_fasta_path:
            ensure_dir(os.path.dirname(out_fasta_path))
            import shutil
            shutil.copyfile(fasta_path, out_fasta_path)
            fasta_path = out_fasta_path

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
        last_err = e

    # 2) Fallback: DFU get_objects (less reliable for assemblies)
    try:
        if dfu is None:
            dfu = DataFileUtil(callback_url)

        obj = dfu.get_objects({"object_refs": [assembly_ref]})["data"][0]["data"]

        contigs: Dict[str, str] = {}
        contig_order: List[str] = []

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
                "or all contigs were filtered out. AssemblyUtil is recommended for Assembly objects."
            )

        # If caller wanted a FASTA written, write it from contigs we recovered
        if out_fasta_path:
            ensure_dir(os.path.dirname(out_fasta_path))
            with open(out_fasta_path, "w") as f:
                for cid in contig_order:
                    f.write(f">{cid}\n")
                    seq = contigs[cid]
                    for i in range(0, len(seq), 60):
                        f.write(seq[i:i+60] + "\n")

        return contig_order, contigs

    except Exception as e:
        raise ValueError(
            "Could not obtain contig sequences from Assembly via AssemblyUtil and DFU fallback. "
            f"AssemblyUtil error was: {repr(last_err)}; DFU error was: {repr(e)}"
        )
