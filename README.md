## Overview

This pipeline processes paired csRNA-seq (capped small RNA) and input (total RNA) libraries to identify and quantify transcription start sites (TSSs) genome-wide. It works with any organism for which you have a genome FASTA, and has been tested with *Arabidopsis thaliana* datasets.

## Installation

### Dependencies

Two environment specs are provided under `workflow/envs/`:

| File | Purpose |
|------|---------|
| `minimal_env.yaml` | Runtime dependencies for running the pipeline with bwa |
| `full_env.yaml` | Everything in minimal plus all aligners, the test suite, and development tools |

Create and activate an environment:

```bash
# For running the pipeline
conda env create -f workflow/envs/minimal_env.yaml
conda activate csRNAseq

# For development and running the full test suite
conda env create -f workflow/envs/full_env.yaml
conda activate csRNAseq_dev
```

At least one dependency will require manual installation, and another must be manually installed on macOS:

- **bfqutils**: See https://github.com/noborilab/bfqutils for installation instructions.
- **HOMER**: If not on Linux, follow the instructions at http://homer.ucsd.edu/homer/introduction/install.html. (To also make use of the environment files, delete HOMER from the yaml file.) After installation, configure the target genome, e.g.:
  ```bash
  perl configureHomer.pl -install hg38    # human
  perl configureHomer.pl -install mm10    # mouse
  perl configureHomer.pl -install tair10  # Arabidopsis
  ```

## Running the Pipeline

```bash
# Minimum 3 cores (some HOMER programs always use 3 threads internally)
snakemake --cores 6 --configfile config/config.yaml

# Dry-run to verify the workflow before executing
snakemake --cores 6 --configfile config/config.yaml -n

```

## Testing

A test suite using a small synthetic genome and synthetic FASTQs lives in `tests/`. It covers schema validation, Snakemake DAG correctness, R-script unit tests, and a full end-to-end run. All committed test data is under 500 KB.

```bash
# Run from the repo root (inside Singularity container or with all tools on PATH)
./tests/run_tests.sh all         # all layers
./tests/run_tests.sh schema      # schema validation only (fastest; no biology tools needed)
./tests/run_tests.sh unit        # R-script unit tests
```

See `tests/README.md` for detailed instructions and how to regenerate fixtures.

## Required config entries

### `sample_table`

A TSV with at least three columns: `sample_name`, `sample_type`, and `read_r1`. Optional columns: `replicate`, `read_r2`, and `input_name`.

- `sample_name`: Base name for this library. csRNA/input pairs share the same `sample_name` by default; replicates of the same condition also share it. Used to group samples during normalization.
- `replicate`: Integer replicate number. Matched csRNA and input libraries must have the same value. Defaults to 1 if the column is absent.
- `sample_type`: `csrna` or `input`.
- `read_r1`: Path(s) to R1 FASTQ files. Comma-separate multiple files for the same sample.
- `read_r2`: Path(s) to R2 FASTQ files for paired-end data. Leave blank for single-end. Mixed single/paired-end datasets are supported by leaving this column empty for the single-end rows.
- `input_name`: On csRNA rows only, the `sample_name` of the input library to pair with. Use this when multiple csRNA conditions share a single input library (leave blank to use the row's own `sample_name`).

```
sample_name	sample_type	read_r1	read_r2	input_name
sample1	csrna	sample1_1.csrna.r1.fq.gz,sample1_2.csrna.r1.fq.gz		
sample1	input	sample1.input.r1.fq.gz	sample1.input.r2.fq.gz	
sample2	csrna	sample2.csrna.r1.fq.gz		sample1
```

### `chrom_sizes`

A two-column (chromosome name, size in bp) TSV **without a header**. A samtools fai index can be used directly. This file controls which chromosomes are kept throughout the pipeline: only chromosomes listed here appear in merged TSSs, the final consensus TSS set, and bigWig tracks. To exclude organellar or other non-nuclear chromosomes, omit them from this file.

### `qc / organelle_chroms`

List of chromosome names to treat as organellar when computing the percentage of nuclear reads (`PctNuclear`). Reads mapping to these chromosomes are subtracted from the total aligned count. Default: `["Pt", "Mt"]` (Arabidopsis). Examples for other organisms:

```yaml
qc:
  organelle_chroms: ["chrM"]          # human / mouse
  organelle_chroms: ["Mt", "Pt"]      # Arabidopsis
  organelle_chroms: []                # skip organelle subtraction
```

### `program / genome_index`

Path to the pre-built genome index for the selected aligner. If the index does not yet exist and `genome_fasta` is provided, the pipeline will build it automatically (stored in `intermediate_dir`). Do not point `genome_index` and `genome_fasta` at the same path; if they match the pipeline assumes the index already exists.

### `program / homer / genome`

HOMER genome name (e.g. `hg38`, `mm10`, `tair10`). Must match a genome configured in your HOMER installation, **or** a file path for a FASTA of your custom genome.

### `filtering / tss_min_cpm` and `filtering / tss_min_samples`

Signal-based filter applied after quantification and before TMM normalization. A TSS is kept in `tss.final.bed` only if its CPM (counts per million library-size reads, without TMM) is ≥ `tss_min_cpm` in at least `tss_min_samples` csRNA samples. Defaults are `tss_min_cpm: 0` and `tss_min_samples: 1`, which keep all TSSs. Example to require ≥ 1 CPM in ≥ 2 samples:

```yaml
filtering:
  tss_min_cpm: 1
  tss_min_samples: 2
```

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
| `collect_consensus_tss` | Build the consensus TSS set (`tss.consensus.bed`) from csRNA samples, requiring detection in at least `tss_min_reps` replicates. TSSs narrower than 150 bp are padded; overlapping TSSs are split. TSSs overlapping miRNA / pre-tRNA loci (if `qc / mirnas` and `qc / trnas` are set) are removed. |
| `quantify_final_tss` | HOMER quantification of the consensus TSSs in all csRNA libraries. |
| `normalize_tss_quantification` | CPM filter (drop TSSs below `filtering / tss_min_cpm` in fewer than `filtering / tss_min_samples` csRNA samples), then TMM normalization (edgeR `TMMwsp`) of the retained counts. Writes the filtered set as `tss.final.bed`. Default CPM threshold is 0, which keeps all TSSs. |
| `quantify_final_in_tss` | HOMER quantification of the filtered `tss.final.bed` in all input libraries. Used by the final QC step. |
| `qc_final_tss` | Compute FRiP, csEnrichment, replicate correlations, and TSS detection rate on the filtered `tss.final.bed`. Writes `qc/qc_final_cs.txt` and `qc/qc_final_in.txt`. |
| `generate_normalized_bw` | Multiply raw bedGraphs by RPM scale factors and export as bigWig. Small RNA regions can optionally be masked. |

## Output Files

All outputs land in `files / output_dir` (default: `results/`):

| File | Description |
|------|-------------|
| `tss.consensus.bed` | Unfiltered consensus TSS set (before CPM filter). |
| `tss.final.bed` | CPM-filtered consensus TSS coordinates (BED6). With default `tss_min_cpm: 0` this equals `tss.consensus.bed`. |
| `tss.final.raw.txt` | Raw tag counts per TSS per csRNA sample (filtered set only). |
| `tss.final.cpm.txt` | TMM-normalized CPM counts. |
| `norm_factors.txt` | edgeR TMM normalization factors and RPM multipliers. |
| `bw/{sample}.rpm.pos.bw` | Forward-strand RPM-normalized bigWig. |
| `bw/{sample}.rpm.neg.bw` | Reverse-strand RPM-normalized bigWig (scores are negative). |
| `qc/qc_cs.txt` | Per-sample QC metrics for csRNA libraries computed on the merged TSS set (FRiP, enrichment, etc.). |
| `qc/qc_in.txt` | Per-sample QC metrics for input libraries computed on the merged TSS set. |
| `qc/qc_final_cs.txt` | Per-sample QC metrics for csRNA libraries re-computed on `tss.final.bed`. FRiP and enrichment here reflect the filtered set used for downstream analysis. |
| `qc/qc_final_in.txt` | Per-sample QC metrics for input libraries re-computed on `tss.final.bed`. |
| `qc/replicate_correlation.txt` | Pairwise Spearman correlation between csRNA replicates of the same `sample_name`. Empty if no `sample_name` has ≥ 2 replicates. |
| `qc/{sample}.trimming.txt` | bfqutils trimming summary. |
| `qc/{sample}.aln.txt` | samtools flagstat alignment summary. |
| `qc/{sample}.aligner.log` | Aligner stderr (bwa/STAR/bowtie2/hisat2). |

## QC Metrics

The `qc_cs.txt` / `qc_in.txt` tables contain:

| Column | Description |
|--------|-------------|
| `csFRiP` | Fraction of nuclear reads falling in csRNA-called TSSs. Should be > 0.9 for a good csRNA library. |
| `sFRiP` | Fraction of nuclear reads in input-called (small RNA) TSSs. |
| `PctNuclear` | Percentage of reads mapping to nuclear chromosomes (i.e. excluding `qc / organelle_chroms`). |
| `csEnrichment` | Ratio of csRNA FRiP to input FRiP — measures how enriched capped initiation signal is relative to the background. |
| `sDepletion` | Inverse ratio of small RNA signal between csRNA and input libraries. |
| `PretRNAPct` | Fraction of input-TSS reads overlapping pre-tRNA loci (requires `qc / trnas`). |
| `PhosEfficiency` | Ratio of pre-tRNA depletion in csRNA vs input — measures 5′-phosphate removal efficiency. |
| `miRNADepletion` | Ratio of miRNA depletion in csRNA vs input — an independent phosphorylation efficiency metric (requires `qc / mirnas`). |
| `StrandBalance` | Fraction of mapped reads on the + strand. Expect ~0.5 in csRNA libraries; large deviations flag adapter contamination, library-prep strand bias, or pile-ups at a few highly expressed loci. Looser bounds in input libraries since small-RNA biology is genuinely strand-skewed. |
| `TSSDetected` | Fraction of merged-set TSSs with ≥ 1 tag in this library. Low values flag undersequenced libraries. |
| `FinalFRiP` | (`qc_final_*` only) FRiP recomputed on `tss.final.bed`. |
| `FinalEnrichment` | (`qc_final_cs.txt` only) csEnrichment recomputed on `tss.final.bed`. |
| `FinalTSSDetected` | (`qc_final_*` only) Fraction of `tss.final.bed` TSSs with ≥ 1 tag in this library. |
| `NConsensusTSS` | (`qc_final_*` only) Number of TSSs in `tss.consensus.bed` (before CPM filter). |
| `NFinalTSS` | (`qc_final_*` only) Number of TSSs retained in `tss.final.bed` after CPM filter. |
| `NFilteredTSS` | (`qc_final_*` only) Number of TSSs dropped by the CPM filter. |
| `MinReplCorrSpearman` | Minimum Spearman correlation between this sample's TSS counts and any other csRNA replicate of the same `sample_name`. `NA` when only one replicate exists. Sharp drops (e.g. < 0.9) flag sample swaps or replicate dropouts. csRNA samples only. |
| `MinReplCorrPearson` | Same as above but Pearson correlation on log1p-transformed counts. |
| `Status` | `Ok` or `FAIL` based on `min_cs_frip` and `min_pct_nuclear` thresholds. |
