# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser

from kb_ConsensusGeneCalling.kb_ConsensusGeneCallingImpl import kb_ConsensusGeneCalling
from kb_ConsensusGeneCalling.kb_ConsensusGeneCallingServer import MethodContext
from kb_ConsensusGeneCalling.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class kb_ConsensusGeneCallingTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_ConsensusGeneCalling'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_ConsensusGeneCalling',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = kb_ConsensusGeneCalling(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_local_prodigal_pipeline_no_services(self):
        # This test bypasses DFU/Report and just validates we can run prodigal + normalize locally.
        module_dir = os.path.dirname(__file__)
        fasta_path = os.path.join(module_dir, "data", "mini.fna")
        self.assertTrue(os.path.exists(fasta_path), f"Missing test FASTA: {fasta_path}")

        from kb_ConsensusGeneCalling.gene_calling_io import stage_input_fasta, ensure_dir
        from kb_ConsensusGeneCalling.callers.prodigal import run_prodigal, parse_fasta_map, parse_prodigal_gff
        from kb_ConsensusGeneCalling.normalize import normalize_feature_objects

        # stage
        scratch = os.path.join("/kb/module/work/tmp", "server_local_mode")
        ensure_dir(scratch)
        staged_fa, contig_order = stage_input_fasta({"input_fasta_path": fasta_path}, scratch)

        # run prodigal
        prod_dir = os.path.join(scratch, "prodigal")
        res = run_prodigal(staged_fa, prod_dir, genetic_code=11, mode="meta")

        raw_feats = parse_prodigal_gff(res["gff"])
        prots = parse_fasta_map(res["faa"])
        dnas = parse_fasta_map(res["ffn"])

        genes = normalize_feature_objects(
            contig_order=contig_order,
            raw_feats=raw_feats,
            proteins=prots,
            dnas=dnas,
            args={
                "coordinate_system": "1-based_closed",
                "min_cds_length_nt": 0,
                "include_partial_genes": True,
                "id_prefix": "ServerLocal"
            }
        )

        self.assertTrue(len(genes) > 0)
        self.assertTrue(genes[0]["ID"].startswith("ServerLocal.peg."))


