"""
Unit test for workflow/scripts/normalize_tss_quantification.R

Run via tests/run_tests.sh or directly:
    tmp=$(mktemp -d)
    snakemake -s tests/unit/normalize.smk --configfile tests/data/config.yaml \
        --config test_outdir=$tmp --cores 1
"""

configfile: "tests/data/config.yaml"

FIXTURES = "tests/data/unit_fixtures"
TEST_OUTDIR = config.get("test_outdir", "tests/unit/tmp_normalize")


rule test_normalize_tss_quantification:
    input:
        quant=FIXTURES + "/tss.consensus.homer.raw.txt",
        bed=FIXTURES + "/tss.consensus.bed",
    output:
        bed=TEST_OUTDIR + "/tss.final.bed",
        raw=TEST_OUTDIR + "/tss.final.raw.txt",
        norm=TEST_OUTDIR + "/tss.final.cpm.txt",
        normfactors=TEST_OUTDIR + "/norm_factors.txt",
    params:
        min_cpm=config["filtering"].get("tss_min_cpm", 0),
        min_samples=config["filtering"].get("tss_min_samples", 1),
    script:
        "../../workflow/scripts/normalize_tss_quantification.R"
