"""
Unit test for workflow/scripts/qc_final_tss.R

Run via tests/run_tests.sh or directly:
    tmp=$(mktemp -d)
    snakemake -s tests/unit/qc_final.smk --configfile tests/data/config.yaml \
        --config test_outdir=$tmp --cores 1
"""

configfile: "tests/data/config.yaml"

FIXTURES = "tests/data/unit_fixtures"
TEST_OUTDIR = config.get("test_outdir", "tests/unit/tmp_qc_final")


CS_IDS = ["condA_csrna1", "condA_csrna2", "condB_csrna1"]
IN_IDS = ["condA_input1", "condA_input2"]


def gate(key, fallback):
    """The same fallback the Snakefile resolves: a final gate left unset uses the
    initial one, so an existing config behaves exactly as it did before the final
    gates existed. An explicit None test, so a gate of 0 is not read as unset."""
    value = config["qc"].get(key)
    return fallback if value is None else value


rule all:
    input:
        TEST_OUTDIR + "/qc_final_cs.txt",
        TEST_OUTDIR + "/qc_final_in.txt",
        TEST_OUTDIR + "/replicate_correlation.txt",


# tss.final.raw.txt must be in the clean format written by normalize_tss_quantification.R:
#   TSS <tab> sample1 <tab> sample2 ...
# We generate it on the fly from the HOMER quant fixture by selecting the keep set
# (all TSSs, since we're not testing the filter here).
# For simplicity we reuse the consensus HOMER quant but expose it via a Snakemake
# rule that converts it to the clean format.

rule make_clean_quant:
    input:
        homer=FIXTURES + "/tss.consensus.homer.raw.txt",
    output:
        clean=TEST_OUTDIR + "/tss.final.raw.txt",
    run:
        import pandas as pd
        df = pd.read_csv(input.homer, sep="\t")
        tag_cols = [c for c in df.columns if " Tag Count" in c]
        counts = df[tag_cols].copy()
        counts.columns = [
            c.split(" Tag Count")[0].rsplit("/", 1)[-1] for c in counts.columns
        ]
        counts.insert(0, "TSS", df.iloc[:, 0])
        counts.to_csv(output.clean, sep="\t", index=False)


rule test_qc_final_tss:
    input:
        quant_cs=TEST_OUTDIR + "/tss.final.raw.txt",
        quant_in=FIXTURES + "/tss.final.in.raw.txt",
        quant_in_cs=FIXTURES + "/quant_in_cs.txt",
        quant_in_in=FIXTURES + "/quant_in_in.txt",
        tss_final=FIXTURES + "/tss.consensus.bed",      # use as proxy for tss.final.bed
        tss_consensus=FIXTURES + "/tss.consensus.bed",
        tss_in=FIXTURES + "/merged_in.bed",
        stats_cs=FIXTURES + "/stats_cs.txt",
        stats_in=FIXTURES + "/stats_in.txt",
        tss_stats=[FIXTURES + "/tss/" + s + ".stats.txt" for s in CS_IDS + IN_IDS],
        repl_cor_initial=FIXTURES + "/repl_cor_initial.txt",
    output:
        qc_cs=TEST_OUTDIR + "/qc_final_cs.txt",
        qc_in=TEST_OUTDIR + "/qc_final_in.txt",
        repl_cor=TEST_OUTDIR + "/replicate_correlation.txt",
    params:
        cs_ids="condA_csrna1 condA_csrna2 condB_csrna1",
        paired_in_ids="condA_input1 condA_input2 condA_input1",
        mirnas=config.get("qc", {}).get("mirnas") or "",
        trnas=config.get("qc", {}).get("trnas") or "",
        min_cs_frip=gate("min_cs_frip_final", config["qc"]["min_cs_frip"]),
        min_pct_nuclear=gate("min_pct_nuclear_final", config["qc"]["min_pct_nuclear"]),
        min_log2_fold=gate(
            "min_log2_fold",
            config["program"]["homer"]["tss"].get("default_log2_fold", 1),
        ),
        tss_dir=FIXTURES + "/tss",
    script:
        "../../workflow/scripts/qc_final_tss.R"
