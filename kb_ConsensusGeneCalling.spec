module kb_ConsensusGeneCalling {

    typedef structure {
        string workspace_name;
        string assembly_ref;
        string output_genome_name;

        string backend;              /* e.g. "consensus" */
        int min_contig_length;
        int genetic_code;
        int min_cds_length_nt;

        int call_cds;                /* 0/1 */
        int call_rnas;               /* 0/1 */
        int include_partial_genes;   /* 0/1 */

        string coordinate_system;    /* "1-based_closed" | "0-based_halfopen" */
    } CallGenesParams;

    typedef structure {
        string report_name;
        string report_ref;
        string genome_ref;
    } CallGenesResult;

    funcdef call_genes(CallGenesParams params) returns (CallGenesResult output);

};
