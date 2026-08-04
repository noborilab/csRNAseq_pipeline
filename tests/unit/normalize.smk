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
        quant_in=FIXTURES + "/tss.consensus.in.homer.raw.txt",
        bed=FIXTURES + "/tss.consensus.bed",
        sizes=FIXTURES + "/tss.consensus.sizes.txt",
    output:
        bed=TEST_OUTDIR + "/tss.final.bed",
        raw=TEST_OUTDIR + "/tss.final.raw.txt",
        norm=TEST_OUTDIR + "/tss.final.cpm.txt",
        normfactors=TEST_OUTDIR + "/norm_factors.txt",
    params:
        min_cpm=config["filtering"].get("tss_min_cpm", 0),
        min_samples=config["filtering"].get("tss_min_samples", 1),
        min_cs_in_ratio=config["filtering"].get("tss_min_cs_in_ratio", 0),
        min_ratio_samples=config["filtering"].get("tss_min_ratio_samples", 1),
        srna_sizes=" ".join(
            str(int(x)) for x in config["filtering"].get("tss_srna_sizes", []) or []
        ),
        max_srna_fraction=config["filtering"].get("tss_max_srna_fraction", 0.5),
        srna_min_reads=config["filtering"].get("tss_srna_min_reads", 30),
        srna_min_samples=config["filtering"].get("tss_srna_min_samples", 1),
        cs_ids="condA_csrna1 condA_csrna2 condB_csrna1",
        paired_in_ids="condA_input1 condA_input2 condA_input1",
    script:
        "../../workflow/scripts/normalize_tss_quantification.R"
