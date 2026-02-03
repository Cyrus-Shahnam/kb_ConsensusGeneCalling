# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import uuid

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil

from kb_ConsensusGeneCalling.gene_calling_io import ensure_dir, fetch_assembly_as_fasta
from kb_ConsensusGeneCalling.callers.prodigal import run_prodigal, parse_fasta_map, parse_prodigal_gff
from kb_ConsensusGeneCalling.callers.trnascan_se import run_trnascan_se, parse_trnascan_tsv
from kb_ConsensusGeneCalling.callers.barrnap import run_barrnap, parse_barrnap_gff

from kb_ConsensusGeneCalling.normalize import normalize_feature_objects
from kb_ConsensusGeneCalling.merge_features import merge_overlapping_features
from kb_ConsensusGeneCalling.write_outputs import write_features_tsv, write_features_json
#END_HEADER


class kb_ConsensusGeneCalling:
    '''
    Module Name:
    kb_ConsensusGeneCalling

    Module Description:
    
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/Cyrus-Shahnam/kb_ConsensusGeneCalling.git"
    GIT_COMMIT_HASH = "c47b3b1aa35a1a828b63f04d465531dbc27b29a8"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.scratch = config.get('scratch', '/kb/module/work/tmp')
        os.makedirs(self.scratch, exist_ok=True)

        self.callback_url = os.environ.get('SDK_CALLBACK_URL')
        if not self.callback_url:
            raise ValueError('SDK_CALLBACK_URL is not set in the environment')

        self.dfu = DataFileUtil(self.callback_url)
        self.report = KBaseReport(self.callback_url)
        #END_CONSTRUCTOR
        pass


    def call_genes(self, ctx, params):
        """
        :param params: instance of type "CallGenesParams" -> structure:
           parameter "workspace_name" of String, parameter "assembly_ref" of
           String, parameter "output_genome_name" of String, parameter
           "backend" of String, parameter "min_contig_length" of Long,
           parameter "genetic_code" of Long, parameter "min_cds_length_nt" of
           Long, parameter "call_cds" of Long, parameter "call_rnas" of Long,
           parameter "include_partial_genes" of Long, parameter
           "coordinate_system" of String
        :returns: instance of type "CallGenesResult" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN call_genes
        workspace_name = params.get("workspace_name")
        assembly_ref = params.get("assembly_ref")
        output_genome_name = params.get("output_genome_name")

        if not workspace_name:
            raise ValueError("workspace_name is required")
        if not assembly_ref:
            raise ValueError("assembly_ref is required")
        if not output_genome_name:
            raise ValueError("output_genome_name is required")

        backend = (params.get("backend") or "prodigal").strip().lower()
        call_cds = bool(params.get("call_cds", True))
        call_rnas = bool(params.get("call_rnas", True))
        min_contig_length = int(params.get("min_contig_length", 200))
        genetic_code = int(params.get("genetic_code", 11))
        min_cds_length_nt = int(params.get("min_cds_length_nt", 90))
        include_partial_genes = bool(params.get("include_partial_genes", False))
        coordinate_system = params.get("coordinate_system", "1-based_closed")

        if coordinate_system not in ("1-based_closed", "0-based_halfopen"):
            raise ValueError(f"Unsupported coordinate_system: {coordinate_system}")

        if backend not in ("prodigal", "prokka"):
            raise ValueError(f"Unsupported backend: {backend}")

        run_id = uuid.uuid4().hex
        run_dir = os.path.join(self.scratch, f"gene_calling_{run_id}")
        ensure_dir(run_dir)

        # ---- Fetch input assembly as FASTA ----
        fasta_path = os.path.join(run_dir, "input.fna")
        contig_order, contigs = fetch_assembly_as_fasta(
            dfu=self.dfu,
            assembly_ref=assembly_ref,
            out_fasta_path=fasta_path,
            min_contig_length=min_contig_length
        )

        notes = []
        all_features = []
        file_links = []

        # ---- CDS: Prodigal ----
        prod_paths = None
        if backend == "prokka":
            notes.append("backend=prokka selected, but current implementation still runs prodigal + RNAs directly (not full prokka).")

        if call_cds:
            prod_out_dir = os.path.join(run_dir, "prodigal")
            prod_paths = run_prodigal(
                fasta_path=fasta_path,
                out_dir=prod_out_dir,
                genetic_code=genetic_code,
                mode="meta"  # safer for fragmented assemblies; avoids 20kb training fail
            )

            raw_feats = parse_prodigal_gff(prod_paths["gff"])
            proteins = parse_fasta_map(prod_paths["faa"])
            dnas = parse_fasta_map(prod_paths["ffn"])

            cds_feats = normalize_feature_objects(
                contig_order=contig_order,
                raw_feats=raw_feats,
                proteins=proteins,
                dnas=dnas,
                args={
                    "coordinate_system": coordinate_system,
                    "min_cds_length_nt": min_cds_length_nt,
                    "include_partial_genes": include_partial_genes,
                    "include_dna_sequence": True,
                    "id_prefix": output_genome_name,
                    "caller": "prodigal",
                }
            )
            all_features.extend(cds_feats)

            file_links.extend([
                {"path": prod_paths["gff"], "name": "genes.gff", "label": "Prodigal GFF"},
                {"path": prod_paths["faa"], "name": "proteins.faa", "label": "Prodigal proteins (FAA)"},
                {"path": prod_paths["ffn"], "name": "genes.ffn", "label": "Prodigal CDS nucleotides (FFN)"},
            ])
        else:
            notes.append("call_cds=false, skipping Prodigal.")

        # ---- RNAs: tRNAscan-SE + barrnap ----
        trna_paths = None
        rrna_paths = None

        if call_rnas:
            # tRNAscan-SE
            trna_out_dir = os.path.join(run_dir, "trnascan_se")
            trna_paths = run_trnascan_se(fasta_path=fasta_path, out_dir=trna_out_dir)
            trna_raw = parse_trnascan_tsv(trna_paths["tsv"])

            trna_feats = normalize_feature_objects(
                contig_order=contig_order,
                raw_feats=trna_raw,
                proteins={},
                dnas={},
                args={
                    "coordinate_system": coordinate_system,
                    "include_dna_sequence": False,
                    "id_prefix": output_genome_name,
                    "caller": "trnascan-se",
                }
            )
            all_features.extend(trna_feats)
            file_links.append({"path": trna_paths["tsv"], "name": "trnascan.tsv", "label": "tRNAscan-SE hits (TSV)"})

            # barrnap
            rrna_out_dir = os.path.join(run_dir, "barrnap")
            rrna_paths = run_barrnap(fasta_path=fasta_path, out_dir=rrna_out_dir)
            rrna_raw = parse_barrnap_gff(rrna_paths["gff"])

            rrna_feats = normalize_feature_objects(
                contig_order=contig_order,
                raw_feats=rrna_raw,
                proteins={},
                dnas={},
                args={
                    "coordinate_system": coordinate_system,
                    "include_dna_sequence": False,
                    "id_prefix": output_genome_name,
                    "caller": "barrnap",
                }
            )
            all_features.extend(rrna_feats)
            file_links.append({"path": rrna_paths["gff"], "name": "barrnap.gff", "label": "barrnap hits (GFF)"})

        else:
            notes.append("call_rnas=false, skipping RNA callers.")

        # ---- Overlap-aware collapsing into regions ----
        regions = merge_overlapping_features(
            features=all_features,
            contig_order=contig_order,
            coordinate_system=coordinate_system,
        )
        
        # ---- Write outputs ----
        features_json_path = os.path.join(run_dir, "features.json")
        features_tsv_path = os.path.join(run_dir, "features.tsv")
        regions_json_path = os.path.join(run_dir, "regions.json")
        regions_tsv_path = os.path.join(run_dir, "regions.tsv")

        write_features_json(all_features, features_json_path)
        write_features_tsv(all_features, features_tsv_path)
        write_features_json(regions, regions_json_path)
        write_features_tsv(regions, regions_tsv_path)

        file_links.extend([
            {"path": features_tsv_path, "name": "features.tsv", "label": "All features (TSV)"},
            {"path": features_json_path, "name": "features.json", "label": "All features (JSON)"},
            {"path": regions_tsv_path, "name": "regions.tsv", "label": "Overlap-collapsed regions (TSV)"},
            {"path": regions_json_path, "name": "regions.json", "label": "Overlap-collapsed regions (JSON)"},
        ])

        # ---- Simple HTML summary ----
        html_path = os.path.join(run_dir, "report.html")
        with open(html_path, "w") as f:
            f.write(f"""<html><body>
<h2>kb_ConsensusGeneCalling</h2>
<ul>
  <li><b>backend</b>: {backend}</li>
  <li><b>assembly_ref</b>: {assembly_ref}</li>
  <li><b>contigs_used</b>: {len(contig_order)} (min_contig_length={min_contig_length})</li>
  <li><b>features_emitted</b>: {len(all_features)}</li>
  <li><b>regions_emitted</b>: {len(regions)}</li>
  <li><b>coordinate_system</b>: {coordinate_system}</li>
</ul>
</body></html>""")

        report_object_name = f"gene_calling_report_{run_id}"
        msg_lines = [
            "Gene calling complete.",
            f"backend={backend}",
            f"features={len(all_features)}",
            f"regions={len(regions)}",
        ]
        if notes:
            msg_lines.append("")
            msg_lines.append("Notes:")
            msg_lines.extend([f"- {n}" for n in notes])

        report_info = self.report.create_extended_report({
            "workspace_name": workspace_name,
            "report_object_name": report_object_name,
            "direct_html_link_index": 0,
            "html_links": [{"path": html_path, "name": "report.html", "label": "Run summary"}],
            "file_links": file_links,
            "message": "\n".join(msg_lines),
        })

        output = {
            "report_name": report_info["name"],
            "report_ref": report_info["ref"],
        }
        #END call_genes

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method call_genes return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
