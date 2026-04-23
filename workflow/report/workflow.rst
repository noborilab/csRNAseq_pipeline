csRNA-seq Analysis Pipeline
============================

This report summarises a run of the csRNA-seq pipeline.

**Alignment program:** ``{{ snakemake.config["program"]["alignment_program"] }}``

**Genome (HOMER):** ``{{ snakemake.config["program"]["homer"]["genome"] }}``

**Minimum replicates for consensus TSSs:** {{ snakemake.config["filtering"]["tss_min_reps"] }}

**QC thresholds:** csRNA FRiP > {{ snakemake.config["qc"]["min_cs_frip"] }}, nuclear reads > {{ snakemake.config["qc"]["min_pct_nuclear"] }}%

Samples processed
-----------------

{% for s in snakemake.config["_samples"] | default([]) %}
- {{ s }}
{% else %}
See the ``qc_cs`` and ``qc_in`` tables for the full sample list.
{% endfor %}

Pipeline stages
---------------

1. Adapter trimming and length truncation (bfqutils)
2. Genome alignment with quality and length filtering
3. HOMER tag directory creation
4. Per-sample TSS detection (``findcsRNATSS.pl``) using paired input library as background
5. Strand-separated raw bedGraph generation
6. Cross-sample QC: FRiP, nuclear read fraction, csRNA enrichment, miRNA/pre-tRNA contamination
7. Consensus TSS set requiring detection in ≥ ``tss_min_reps`` replicates
8. TMM normalisation (edgeR ``TMMwsp``)
9. RPM-normalised bigWig track generation
