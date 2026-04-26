## Overview

This pipeline processes paired csRNA-seq (capped small RNA) and input (total RNA) libraries to identify and quantify transcription start sites (TSSs) genome-wide. It was developed for *Arabidopsis thaliana* but is usable with any organism whose genome is configured in HOMER.

## Installation

### Dependencies

Install the conda environment:

```bash
conda env create -f workflow/envs/env.yaml
conda activate csRNAseq_snakemake
```

Two tools **must be installed manually** before running the pipeline — they are not on Bioconda:

- **HOMER** — follow the instructions at http://homer.ucsd.edu/homer/introduction/install.html. After installation, configure the target genome, e.g.:
  ```bash
  perl configureHomer.pl -install tair10
  ```
- **bfqutils** — see the project repository for installation instructions.

## Running the Pipeline

```bash
# Minimum 3 cores (some HOMER programs always use 3 threads internally)
snakemake --cores 6 --configfile config/config.yaml

# Dry-run to verify the workflow before executing
snakemake --cores 6 --configfile config/config.yaml -n

```

## Required config entries

### `sample_table`

A TSV with at least three columns: `sample_name`, `sample_type`, and `read_r1`. Optional columns: `replicate` and `read_r2`.

- `sample_name`: Base name shared by a matched csRNA/input pair (and replicates of the same condition). Used to group samples during normalization.
- `replicate`: Integer replicate number. Matched csRNA and input libraries must have the same value. Defaults to 1 if the column is absent.
- `sample_type`: `csrna` or `input`.
- `read_r1`: Path(s) to R1 FASTQ files. Comma-separate multiple files for the same sample.
- `read_r2`: Path(s) to R2 FASTQ files for paired-end data. Leave blank for single-end. Mixed single/paired-end datasets are supported by leaving this column empty for the single-end rows.

```
sample_name	sample_type	read_r1	read_r2
sample1	csrna	sample1_1.csrna.r1.fq.gz,sample1_2.csrna.r1.fq.gz	
sample1	input	sample1.input.r1.fq.gz	sample1.input.r2.fq.gz
```

### `chrom_sizes`

A two-column (chromosome name, size in bp) TSV **without a header**. A samtools fai index can be used directly. This file controls which chromosomes are kept throughout the pipeline: only chromosomes listed here appear in merged TSSs, the final consensus TSS set, and bigWig tracks. To exclude organellar chromosomes (e.g. `Mt`, `Pt` in Arabidopsis; `chrM` in human), simply omit them from this file.

### `program / genome_index`

Path to the pre-built genome index for the selected aligner. If the index does not yet exist and `genome_fasta` is provided, the pipeline will build it automatically (stored in `intermediate_dir`). Do not point `genome_index` and `genome_fasta` at the same path; if they match the pipeline assumes the index already exists.

### `program / homer / genome`

HOMER genome name (e.g. `tair10`, `hg38`). Must match a genome configured in your HOMER installation.

## Workflow Steps

| Rule | Description |
|------|-------------|
| `build_index` | (Conditional) Build a genome index for the selected aligner from `genome_fasta` if no index is found at `genome_index`. |
| `trim` | Adapter trimming and length truncation (bfqutils). Stats written to `qc/`. |
| `align` | Genome alignment (bwa-aln / bwa-mem / STAR / bowtie2 / hisat2), MAPQ and length filtering, BAM sorting and indexing. |
| `make_tagdir` | HOMER tag directory from BAM. |
| `find_tss_initial` | Per-sample TSS calling (HOMER `findcsRNATSS.pl`) using the paired input library as background. Run on both csRNA and input samples. |
| `make_raw_bedgraph` | Strand-separated raw bedGraph files (HOMER `makeUCSCfile`). |
| `gather_stats` | Extract total reads, organellar reads, and tag frequencies from HOMER tagInfo files. |
| `merge_initial_tss` | Strand-aware merge of all per-sample TSS BEDs. Chromosomes not in `chrom_sizes` are excluded. |
| `quantify_initial_{cs,in}_tss` | HOMER `annotatePeaks.pl` quantification of merged TSS sets in all libraries. |
| `qc_initial_tss` | Calculate FRiP, nuclear read %, csRNA enrichment, miRNA depletion, and phosphorylation efficiency. Samples with FRiP < `min_cs_frip` or nuclear % < `min_pct_nuclear` are flagged as FAIL. |
| `collect_consensus_tss` | Build the final consensus TSS set from csRNA samples, requiring detection in at least `tss_min_reps` replicates. TSSs narrower than 150 bp are padded; overlapping TSSs are split. TSSs overlapping miRNA / pre-tRNA loci (if `qc / mirnas` and `qc / trnas` are set) are removed. |
| `quantify_final_tss` | HOMER quantification of the consensus TSSs in all csRNA libraries. |
| `normalize_tss_quantification` | TMM normalization (edgeR `TMMwsp`) of raw tag counts. |
| `generate_normalized_bw` | Multiply raw bedGraphs by RPM scale factors and export as bigWig. Small RNA regions can optionally be masked. |

## Output Files

All outputs land in `files / output_dir` (default: `results/`):

| File | Description |
|------|-------------|
| `tss.final.bed` | Consensus TSS coordinates (BED6). |
| `tss.final.raw.txt` | Raw tag counts per TSS per csRNA sample. |
| `tss.final.cpm.txt` | TMM-normalized CPM counts. |
| `norm_factors.txt` | edgeR TMM normalization factors and RPM multipliers. |
| `bw/{sample}.rpm.pos.bw` | Forward-strand RPM-normalized bigWig. |
| `bw/{sample}.rpm.neg.bw` | Reverse-strand RPM-normalized bigWig (scores are negative). |
| `qc/qc_cs.txt` | Per-sample QC metrics for csRNA libraries (FRiP, enrichment, etc.). |
| `qc/qc_in.txt` | Per-sample QC metrics for input libraries. |
| `qc/{sample}.trimming.txt` | bfqutils trimming summary. |
| `qc/{sample}.aln.txt` | samtools flagstat alignment summary. |
| `qc/{sample}.aligner.log` | Aligner stderr (bwa/STAR/bowtie2/hisat2). |

## QC Metrics

The `qc_cs.txt` / `qc_in.txt` tables contain:

| Column | Description |
|--------|-------------|
| `csFRiP` | Fraction of nuclear reads falling in csRNA-called TSSs. Should be > 0.9 for a good csRNA library. |
| `sFRiP` | Fraction of nuclear reads in input-called (small RNA) TSSs. |
| `PctNuclear` | Percentage of reads mapping to nuclear chromosomes. |
| `csEnrichment` | Ratio of csRNA FRiP to input FRiP — measures how enriched capped initiation signal is relative to the background. |
| `sDepletion` | Inverse ratio of small RNA signal between csRNA and input libraries. |
| `PretRNAPct` | Fraction of input-TSS reads overlapping pre-tRNA loci (requires `qc / trnas`). |
| `PhosEfficiency` | Ratio of pre-tRNA depletion in csRNA vs input — measures 5′-phosphate removal efficiency. |
| `miRNADepletion` | Ratio of miRNA depletion in csRNA vs input — an independent phosphorylation efficiency metric (requires `qc / mirnas`). |
| `Status` | `Ok` or `FAIL` based on `min_cs_frip` and `min_pct_nuclear` thresholds. |
