# -*- coding: utf-8 -*-
import os
import unittest

from kb_ConsensusGeneCalling.callers.trnascan_se import parse_trnascan_tsv
from kb_ConsensusGeneCalling.callers.barrnap import parse_barrnap_gff


class TestRNAParsers(unittest.TestCase):

    def test_parse_trnascan_tsv(self):
        module_dir = os.path.dirname(__file__)
        tsv_path = os.path.join(module_dir, "data", "trnascan_example.tsv")
        feats = parse_trnascan_tsv(tsv_path)

        self.assertEqual(len(feats), 2)

        f0 = feats[0]
        self.assertEqual(f0["contig"], "contigA")
        self.assertEqual(f0["start"], 10)
        self.assertEqual(f0["stop"], 50)
        self.assertEqual(f0["strand"], "+")
        self.assertEqual(f0["type"], "tRNA")
        self.assertTrue(f0["id"].startswith("contigA.tRNA."))

        f1 = feats[1]
        self.assertEqual(f1["contig"], "contigB")
        self.assertEqual(f1["start"], 80)
        self.assertEqual(f1["stop"], 120)
        self.assertEqual(f1["strand"], "-")
        self.assertEqual(f1["type"], "tRNA")

    def test_parse_barrnap_gff(self):
        module_dir = os.path.dirname(__file__)
        gff_path = os.path.join(module_dir, "data", "barrnap_example.gff")
        feats = parse_barrnap_gff(gff_path)

        self.assertEqual(len(feats), 2)

        f0 = feats[0]
        self.assertEqual(f0["contig"], "contigA")
        self.assertEqual(f0["start"], 5)
        self.assertEqual(f0["stop"], 120)
        self.assertEqual(f0["strand"], "+")
        self.assertEqual(f0["type"], "rRNA")
        self.assertEqual(f0["id"], "rrna1")

        f1 = feats[1]
        self.assertEqual(f1["contig"], "contigB")
        self.assertEqual(f1["strand"], "-")
        self.assertEqual(f1["type"], "rRNA")
        self.assertEqual(f1["id"], "rrna2")
