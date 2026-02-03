import os
import subprocess
import gzip
import shutil
from typing import Dict, List, Tuple

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil


def parse_fasta_ids(fasta_path: str) -> List[str]:
    """
    Return a list of FASTA record IDs in file order.
    ID is the header up to first whitespace (matches _read_fasta_records).
    """
    return [cid for cid, _seq in _read_fasta_records(fasta_path)]


def _copy_or_decompress_to(src_path: str, dst_path: str) -> str:
    """
    Copy src_path to dst_path. If src_path ends with .gz, decompress while copying.
    Returns dst_path.
    """
    ensure_dir(os.path.dirname(dst_path) or ".")
    if src_path.endswith(".gz"):
        with gzip.open(src_path, "rt") as fin, open(dst_path, "w") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copyfile(src_path, dst_path)
    return dst_path


def stage_input_fasta(
    *,
    assembly_ref: str = None,
    fasta_path: str = None,
    out_fasta_path: str,
    min_contig_length: int = 0,
    dfu=None,
    ctx=None,
    scratch: str = None,
) -> Tuple[List[str], Dict[str, str], str]:
    """
    Stage an input FASTA for downstream callers.

    Exactly one of (assembly_ref, fasta_path) must be provided.

    - If assembly_ref is provided: uses fetch_assembly_as_fasta to write a FASTA.
    - If fasta_path is provided: copies (or gunzips) into out_fasta_path.
    - Always loads contigs into memory and returns (contig_order, contigs, out_fasta_path).

    Returns:
      contig_order: list[str]
      contigs: dict[str,str]
      staged_fasta_path: str
    """
    if bool(assembly_ref) == bool(fasta_path):
        raise ValueError("Provide exactly one of assembly_ref or fasta_path")

    ensure_dir(os.path.dirname(out_fasta_path) or ".")

    if assembly_ref:
        contig_order, contigs = fetch_assembly_as_fasta(
            ctx=ctx,
            assembly_ref=assembly_ref,
            scratch=scratch,
            min_contig_length=min_contig_length,
            dfu=dfu,
            out_fasta_path=out_fasta_path,
        )
        return contig_order, contigs, out_fasta_path

    # local FASTA input
    _copy_or_decompress_to(fasta_path, out_fasta_path)
    contig_order, contigs = _load_contigs_from_fasta(out_fasta_path, min_contig_length=min_contig_length)
    return contig_order, contigs, out_fasta_path

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


def _unwrap_and_validate_return(contig_order, contigs):
    """
    Some call sites / refactors accidentally produce contig_order as (order, contigs).
    This unwraps that case and validates we return (list[str], dict[str,str]).
    """
    # Unwrap if contig_order is actually a packed (order, contigs)
    if (
        isinstance(contig_order, tuple)
        and len(contig_order) == 2
        and isinstance(contig_order[0], list)
        and isinstance(contig_order[1], dict)
        and contigs is not None
        and isinstance(contigs, dict)
        and len(contigs) == 0
    ):
        # In this pattern, contigs arg is often empty while packed tuple holds the real dict
        contig_order, contigs = contig_order

    # Validate types
    if not isinstance(contig_order, list):
        raise TypeError(f"contig_order must be list[str], got {type(contig_order)}: {contig_order}")

    if not contig_order:
        raise ValueError("contig_order is empty (no contigs after filtering)")

    bad = [x for x in contig_order[:10] if not isinstance(x, str)]
    if bad:
        raise TypeError(f"contig_order must be list[str]; found {type(bad[0])} value={bad[0]}")

    if not isinstance(contigs, dict):
        raise TypeError(f"contigs must be dict[str,str], got {type(contigs)}")

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
    if not assembly_ref:
        raise ValueError("assembly_ref is required")

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

        contig_order, contigs = _unwrap_and_validate_return(contig_order, contigs)
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
                        f.write(seq[i:i + 60] + "\n")

        contig_order, contigs = _unwrap_and_validate_return(contig_order, contigs)
        return contig_order, contigs

    except Exception as e:
        raise ValueError(
            "Could not obtain contig sequences from Assembly via AssemblyUtil and DFU fallback. "
            f"AssemblyUtil error was: {repr(last_err)}; DFU error was: {repr(e)}"
        )
