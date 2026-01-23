# -*- coding: utf-8 -*-
import os
import unittest

from kb_ConsensusGeneCalling.callers.prodigal import run_prodigal, parse_prodigal_gff, parse_fasta_map
from kb_ConsensusGeneCalling.gene_calling_io import ensure_dir, parse_fasta_ids
from kb_ConsensusGeneCalling.normalize import normalize_feature_objects


class TestNormalize(unittest.TestCase):

    def test_normalize_coordinate_systems_and_filters(self):
        module_dir = os.path.dirname(__file__)
        fasta_path = os.path.join(module_dir, "data", "mini.fna")
        self.assertTrue(os.path.exists(fasta_path), f"Missing test FASTA: {fasta_path}")

        contig_order = parse_fasta_ids(fasta_path)
        self.assertTrue(len(contig_order) > 0, "Expected contigs in mini.fna")

        out_dir = os.path.join("/kb/module/work/tmp", "ut_normalize")
        ensure_dir(out_dir)

        res = run_prodigal(fasta_path=fasta_path, out_dir=out_dir, genetic_code=11, mode="meta")
        raw_feats = parse_prodigal_gff(res["gff"])
        prots = parse_fasta_map(res["faa"])
        dnas = parse_fasta_map(res["ffn"])

        # 1) 1-based closed
        args_1based = {
            "coordinate_system": "1-based_closed",
            "include_dna_sequence": True,
            "min_cds_length_nt": 0,
            "include_partial_genes": True,
            "id_prefix": "TestGenome"
        }
        genes_1 = normalize_feature_objects(contig_order, raw_feats, prots, dnas, args_1based)
        self.assertTrue(len(genes_1) > 0, "Expected normalized CDS features (1-based)")

        g0 = genes_1[0]
        self.assertTrue(g0["start"] >= 1)
        self.assertTrue(g0["stop"] >= 1)
        self.assertIn(g0["direction"], ["+", "-"])
        self.assertTrue(g0["ID"].startswith("TestGenome.peg."))

        # 2) 0-based half-open: start should be 1 less than 1-based for same feature
        args_0based = dict(args_1based)
        args_0based["coordinate_system"] = "0-based_halfopen"
        genes_0 = normalize_feature_objects(contig_order, raw_feats, prots, dnas, args_0based)
        self.assertEqual(len(genes_0), len(genes_1), "Expected same number of genes under coordinate transform")

        # Compare first gene (ordering deterministic)
        self.assertEqual(genes_0[0]["stop"], genes_1[0]["stop"])
        self.assertEqual(genes_0[0]["start"], genes_1[0]["start"] - 1)

        # 3) min length filter should reduce or equal count
        args_minlen = dict(args_1based)
        args_minlen["min_cds_length_nt"] = 999999  # absurdly high -> should remove all
        genes_min = normalize_feature_objects(contig_order, raw_feats, prots, dnas, args_minlen)
        self.assertEqual(len(genes_min), 0, "Expected zero genes after huge min length filter")
