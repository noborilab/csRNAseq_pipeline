"""
Unit test for workflow/scripts/tss_size_composition.R

Run via tests/run_tests.sh or directly:
    tmp=$(mktemp -d)
    snakemake -s tests/unit/size_composition.smk --configfile tests/data/config.yaml \
        --config test_outdir=$tmp --cores 1
"""

configfile: "tests/data/config.yaml"

FIXTURES = "tests/data/unit_fixtures"
TEST_OUTDIR = config.get("test_outdir", "tests/unit/tmp_size_composition")


rule test_tss_size_composition:
    input:
        bed=FIXTURES + "/tss.consensus.bed",
        td=[
            FIXTURES + "/tagdir/condA_csrna1",
            FIXTURES + "/tagdir/condA_csrna2",
            FIXTURES + "/tagdir/condB_csrna1",
        ],
    output:
        TEST_OUTDIR + "/tss.consensus.sizes.txt",
    params:
        srna_sizes=" ".join(
            str(int(x)) for x in config["filtering"].get("tss_srna_sizes", []) or []
        ),
    script:
        "../../workflow/scripts/tss_size_composition.R"
