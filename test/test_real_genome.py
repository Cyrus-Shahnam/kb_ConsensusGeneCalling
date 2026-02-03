# -*- coding: utf-8 -*-
import os, gzip, shutil
from pathlib import Path
import unittest

from kb_ConsensusGeneCalling.callers.prodigal import run_prodigal
from kb_ConsensusGeneCalling.callers.prodigal import parse_prodigal_gff

def gunzip(gz_path: str, out_path: str) -> str:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with gzip.open(gz_path, "rb") as fin, open(out_path, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    return out_path

class TestRealGenome(unittest.TestCase):
    def test_myco_g37_single_mode(self):
        test_dir = Path(__file__).resolve().parent
        fasta = test_dir / "data" / "myco_g37.fna"

#        fasta = gunzip(str(gz), "/kb/module/work/tmp/ut_myco_g37/myco_g37.fna")
        out_dir = "/kb/module/work/tmp/ut_myco_g37/out"

        res = run_prodigal(fasta_path=fasta, out_dir=out_dir, genetic_code=11, mode="single")
        genes = parse_prodigal_gff(res["gff"])

        self.assertGreater(len(genes), 100)   # stable sanity check
