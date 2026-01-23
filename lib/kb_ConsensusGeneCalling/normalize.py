# -*- coding: utf-8 -*-
from __future__ import annotations


def _cds_len_nt(start_1based: int, stop_1based: int) -> int:
    return abs(int(stop_1based) - int(start_1based)) + 1


def normalize_feature_objects(
    contig_order: list[str],
    raw_feats: list[dict],
    proteins: dict[str, str],
    dnas: dict[str, str],
    args: dict
) -> list[dict]:
    """
    Produces standardized feature schema.

    Coordinates:
      - raw_feats are 1-based closed (common for GFF callers like prodigal/barrnap)
      - output controlled by args["coordinate_system"]:
          * 1-based_closed: start/stop inclusive
          * 0-based_halfopen: [start-1, stop) (stop is exclusive)
    """
    coordinate_system = args.get("coordinate_system", "1-based_closed")
    include_dna_sequence = bool(args.get("include_dna_sequence", True))

    min_cds_length_nt = int(args.get("min_cds_length_nt", 0))
    include_partial_genes = bool(args.get("include_partial_genes", False))

    prefix = args.get("id_prefix") or args.get("genome_id") or "genome"
    caller = args.get("caller", "unknown")

    contig_rank = {cid: i for i, cid in enumerate(contig_order)}
    raw_feats_sorted = sorted(
        raw_feats,
        key=lambda r: (contig_rank.get(r["contig"], 10**9), int(r["start"]), int(r["stop"]))
    )

    out: list[dict] = []
    counters = {"CDS": 0, "tRNA": 0, "rRNA": 0, "ncRNA": 0}

    for rf in raw_feats_sorted:
        contig = rf["contig"]
        start = int(rf["start"])
        stop = int(rf["stop"])
        strand = rf["strand"]

        ftype = rf.get("type", "CDS")
        feat_id = rf.get("prodigal_id") or rf.get("id") or ""

        # CDS length filter only
        if ftype == "CDS":
            length_nt = _cds_len_nt(start, stop)
            if length_nt < min_cds_length_nt:
                continue

        prot_seq = proteins.get(feat_id, "") if ftype == "CDS" else ""
        dna_seq = dnas.get(feat_id, "") if include_dna_sequence else ""

        # Partial filtering (CDS only)
        if ftype == "CDS":
            is_partial = False
            if rf.get("partial") is not None:
                v = str(rf["partial"]).strip()
                # Prodigal meta often uses partial=11,00 etc. Treat 00 as "not partial".
                is_partial = (v != "00")
            else:
                if not prot_seq:
                    is_partial = True

            if is_partial and not include_partial_genes:
                continue

        # Coordinate transform
        if coordinate_system == "1-based_closed":
            out_start, out_stop = start, stop
        elif coordinate_system == "0-based_halfopen":
            out_start, out_stop = start - 1, stop
        else:
            raise ValueError(f"Unknown coordinate_system: {coordinate_system}")

        counters.setdefault(ftype, 0)
        counters[ftype] += 1
        n = counters[ftype]

        if ftype == "CDS":
            out_id = f"{prefix}.peg.{n}"
        else:
            out_id = f"{prefix}.{ftype}.{n}"

        feature_obj = {
            "caller": caller,
            "ID": out_id,
            "contig": contig,
            "start": out_start,
            "stop": out_stop,
            "direction": "+" if strand == "+" else "-",
            "DNA sequence": dna_seq,
            "protein sequence": prot_seq if ftype == "CDS" else "",
            "type": ftype
        }

        out.append(feature_obj)

    return out
