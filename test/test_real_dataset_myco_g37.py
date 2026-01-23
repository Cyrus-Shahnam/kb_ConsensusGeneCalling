import os
import gzip
import shutil
import uuid
import unittest

from kb_ConsensusGeneCalling.kb_ConsensusGeneCallingImpl import kb_ConsensusGeneCalling
from installed_clients.WorkspaceClient import Workspace
from installed_clients.DataFileUtilClient import DataFileUtil


class RealDatasetMycoG37SmokeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module_dir = "/kb/module"
        cls.data_dir = os.path.join(cls.module_dir, "test", "data")
        cls.scratch = os.path.join(cls.module_dir, "workdir")
        os.makedirs(cls.scratch, exist_ok=True)

        # Impl + config (deploy.cfg provides workspace-url; test harness sets SDK_CALLBACK_URL + token)
        cls.impl = kb_ConsensusGeneCalling({"scratch": cls.scratch})
        cls.config = cls.impl.config

        cls.ws = Workspace(cls.config["workspace-url"])
        cls.dfu = DataFileUtil(os.environ["SDK_CALLBACK_URL"])

        cls.ws_name = "kb_congenecalling_test_" + str(uuid.uuid4())
        cls.ws.create_workspace({"workspace": cls.ws_name})

    @classmethod
    def tearDownClass(cls):
        # Keep things clean
        try:
            cls.ws.delete_workspace({"workspace": cls.ws_name})
        except Exception:
            pass

    def _gunzip(self, gz_path, out_path):
        with gzip.open(gz_path, "rb") as fin, open(out_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)

    def test_call_genes_on_myco_g37(self):
        gz_fna = os.path.join(self.data_dir, "myco_g37.fna.gz")
        self.assertTrue(os.path.exists(gz_fna), f"Missing dataset: {gz_fna}")

        fna = os.path.join(self.scratch, "myco_g37.fna")
        self._gunzip(gz_fna, fna)
        self.assertGreater(os.path.getsize(fna), 0)

        # Import FASTA -> Assembly (returns an assembly_ref in the workspace)
        assembly_name = "myco_g37_assembly_" + str(uuid.uuid4())[:8]
        imp = self.dfu.fasta_to_assembly({
            "file": {"path": fna},
            "workspace_name": self.ws_name,
            "assembly_name": assembly_name
        })
        assembly_ref = imp["assembly_ref"]

        # Call the KBase method using EXACT param names from spec.json mapping
        params = {
            "workspace_name": self.ws_name,
            "assembly_ref": assembly_ref,
            "output_genome_name": "GeneCalling_" + str(uuid.uuid4())[:8],
            "backend": "prokka",          # "prokka" runs prodigal + trnascan + barrnap
            "call_cds": True,
            "call_rnas": True,
            "min_contig_length": 200,
            "genetic_code": 11,
            "min_cds_length_nt": 90,
            "include_partial_genes": False,
            "coordinate_system": "1-based_closed",
        }

        ctx = None
        result = self.impl.call_genes(ctx, params)

        # Your spec output mapping expects report_name and report_ref
        self.assertIsInstance(result, list)
        self.assertTrue(result and isinstance(result[0], dict))
        self.assertIn("report_name", result[0])
        self.assertIn("report_ref", result[0])

        print("REPORT:", result[0]["report_name"], result[0]["report_ref"])
        print("WORKSPACE:", self.ws_name)
        print("ASSEMBLY_REF:", assembly_ref)
