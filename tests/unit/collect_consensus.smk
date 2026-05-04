"""
Unit test for workflow/scripts/collect_consensus_tss.R

Run via tests/run_tests.sh or directly:
    tmp=$(mktemp -d)
    snakemake -s tests/unit/collect_consensus.smk --configfile tests/data/config.yaml \
        --config test_outdir=$tmp --cores 1
"""

configfile: "tests/data/config.yaml"

FIXTURES = "tests/data/unit_fixtures"
TEST_OUTDIR = config.get("test_outdir", "tests/unit/tmp_collect")


rule test_collect_consensus_tss:
    input:
        tss_cs=[
            FIXTURES + "/per_sample_tss/condA_csrna1.tss.bed",
            FIXTURES + "/per_sample_tss/condA_csrna2.tss.bed",
            FIXTURES + "/per_sample_tss/condB_csrna1.tss.bed",
        ],
        tss_in=FIXTURES + "/merged_in.bed",
    output:
        TEST_OUTDIR + "/tss.consensus.bed",
    params:
        min_reps=config["filtering"]["tss_min_reps"],
        ids="condA_csrna1 condA_csrna2 condB_csrna1 condA_input1 condA_input2",
        sample_names="condA condA condB condA condA",
        sample_types="csrna csrna csrna input input",
        replicates="1 2 1 1 2",
        bed_dir=FIXTURES + "/per_sample_tss",
        bg_dir=TEST_OUTDIR + "/bedgraph",
        chrom_sizes="tests/data/chrom.sizes",
        mirnas="tests/data/mirnas.bed",
        trnas="tests/data/trnas.bed",
    script:
        "../../workflow/scripts/collect_consensus_tss.R"
