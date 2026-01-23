# -*- coding: utf-8 -*-
import os
import unittest

from kb_ConsensusGeneCalling.callers.prodigal import run_prodigal, parse_prodigal_gff, parse_fasta_map
from kb_ConsensusGeneCalling.gene_calling_io import ensure_dir


class TestProdigalParse(unittest.TestCase):

    def test_run_and_parse_prodigal(self):
        module_dir = os.path.dirname(__file__)
        fasta_path = os.path.join(module_dir, "data", "mini.fna")
        self.assertTrue(os.path.exists(fasta_path), f"Missing test FASTA: {fasta_path}")

        out_dir = os.path.join("/kb/module/work/tmp", "ut_prodigal_parse")
        ensure_dir(out_dir)

        res = run_prodigal(
            fasta_path=fasta_path,
            out_dir=out_dir,
            genetic_code=11,
            mode="meta"
        )

        self.assertTrue(os.path.exists(res["gff"]), "Missing genes.gff")
        self.assertTrue(os.path.exists(res["faa"]), "Missing proteins.faa")
        self.assertTrue(os.path.exists(res["ffn"]), "Missing genes.ffn")

        feats = parse_prodigal_gff(res["gff"])
        self.assertTrue(len(feats) > 0, "Expected at least 1 CDS feature in parsed GFF")

        # sanity: parsed coordinates are ints and 1-based inclusive
        f0 = feats[0]
        self.assertIn("contig", f0)
        self.assertIn("start", f0)
        self.assertIn("stop", f0)
        self.assertIn("strand", f0)
        self.assertIn("prodigal_id", f0)
        self.assertIsInstance(f0["start"], int)
        self.assertIsInstance(f0["stop"], int)
        self.assertTrue(f0["start"] >= 1)
        self.assertTrue(f0["stop"] >= 1)

        # sanity: FAA/FFN parse yields some sequences
        prots = parse_fasta_map(res["faa"])
        dnas = parse_fasta_map(res["ffn"])
        self.assertTrue(len(prots) > 0, "Expected proteins in FAA")
        self.assertTrue(len(dnas) > 0, "Expected genes in FFN")
