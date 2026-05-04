"""
Unit test for workflow/scripts/qc_initial_tss.R

Exercises the shared-input pairing (condB_csrna1 → condA_input1 via input_name).

Run via tests/run_tests.sh or directly:
    tmp=$(mktemp -d)
    snakemake -s tests/unit/qc_initial.smk --configfile tests/data/config.yaml \
        --config test_outdir=$tmp --cores 1
"""

configfile: "tests/data/config.yaml"

FIXTURES = "tests/data/unit_fixtures"
TEST_OUTDIR = config.get("test_outdir", "tests/unit/tmp_qc_initial")


rule test_qc_initial_tss:
    input:
        quant_cs_cs=FIXTURES + "/tss.consensus.homer.raw.txt",
        quant_cs_in=FIXTURES + "/quant_cs_in.txt",
        quant_in_cs=FIXTURES + "/quant_in_cs.txt",
        quant_in_in=FIXTURES + "/quant_in_in.txt",
        stats_cs=FIXTURES + "/stats_cs.txt",
        stats_in=FIXTURES + "/stats_in.txt",
        tss_cs=FIXTURES + "/merged_cs.bed",
        tss_in=FIXTURES + "/merged_in.bed",
    output:
        qc_cs=TEST_OUTDIR + "/qc_cs.txt",
        qc_in=TEST_OUTDIR + "/qc_in.txt",
        repl_cor=TEST_OUTDIR + "/replicate_correlation.txt",
    params:
        mirnas=config.get("qc", {}).get("mirnas") or "",
        trnas=config.get("qc", {}).get("trnas") or "",
        # condA_csrna1→condA_input1, condA_csrna2→condA_input2,
        # condB_csrna1→condA_input1 (shared via input_name)
        cs_ids="condA_csrna1 condA_csrna2 condB_csrna1",
        paired_in_ids="condA_input1 condA_input2 condA_input1",
    script:
        "../../workflow/scripts/qc_initial_tss.R"
