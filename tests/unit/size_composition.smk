"""
Unit tests for workflow/scripts/tss_size_composition.R (per library) and
workflow/scripts/merge_size_composition.R (the join).

Run via tests/run_tests.sh or directly:
    tmp=$(mktemp -d)
    snakemake -s tests/unit/size_composition.smk --configfile tests/data/config.yaml \
        --config test_outdir=$tmp --cores 1
"""

configfile: "tests/data/config.yaml"

FIXTURES = "tests/data/unit_fixtures"
TEST_OUTDIR = config.get("test_outdir", "tests/unit/tmp_size_composition")
LIBS = ["condA_csrna1", "condA_csrna2", "condB_csrna1"]


rule all_size_composition:
    input:
        TEST_OUTDIR + "/tss.consensus.sizes.txt",


rule test_tss_size_composition:
    input:
        bed=FIXTURES + "/tss.consensus.bed",
        td=FIXTURES + "/tagdir/{lib}",
    output:
        TEST_OUTDIR + "/tss_sizes/{lib}.sizes.txt",
    params:
        srna_sizes=" ".join(
            str(int(x)) for x in config["filtering"].get("tss_srna_sizes", []) or []
        ),
        top_sizes_max=config["filtering"].get("tss_top_sizes_max", 5),
        max_top_sizes_fraction=config["filtering"].get("tss_max_top_sizes_fraction", 1),
    script:
        "../../workflow/scripts/tss_size_composition.R"


rule test_merge_size_composition:
    input:
        bed=FIXTURES + "/tss.consensus.bed",
        sizes=[TEST_OUTDIR + "/tss_sizes/" + lib + ".sizes.txt" for lib in LIBS],
    output:
        TEST_OUTDIR + "/tss.consensus.sizes.txt",
    script:
        "../../workflow/scripts/merge_size_composition.R"
