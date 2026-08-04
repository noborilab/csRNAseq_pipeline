## Overview

This pipeline takes paired csRNA-seq (capped small RNA) and input (total RNA) libraries and uses them to find and quantify transcription start sites (TSSs) across the genome. In principle it should work with any organism for which you have a genome FASTA, though so far it has really only been tested on *Arabidopsis thaliana* data.

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

A few of the dependencies are not on conda (or not for every platform), so they need to be installed by hand:

- **bfqutils**: see https://github.com/noborilab/bfqutils for installation instructions.
- **bam2td**: see https://github.com/noborilab/bam2td for installation instructions.
- **HOMER**: if you are not on Linux, follow the instructions at http://homer.ucsd.edu/homer/introduction/install.html. (If you would still like to use the environment files in that case, just delete the HOMER line from the yaml first.) Once it is installed, configure the genome you are working with, e.g.:
  ```bash
  perl configureHomer.pl -install hg38    # human
  perl configureHomer.pl -install mm10    # mouse
  perl configureHomer.pl -install tair10  # Arabidopsis
  ```

## Running the Pipeline

With the environment set up and your config and sample table filled in, running the pipeline is just an ordinary Snakemake invocation. You will want at least three cores, since a couple of the HOMER programs always run three threads internally.

```bash
# Minimum 3 cores (some HOMER programs always use 3 threads internally)
snakemake --cores 6 --configfile config/config.yaml

# A dry run first, to check the workflow before actually executing it
snakemake --cores 6 --configfile config/config.yaml -n

```

## Testing

There is a test suite in `tests/`, built around a small synthetic genome and a set of synthetic FASTQs. It covers schema validation, the correctness of the Snakemake DAG, the R scripts on their own, and a full end-to-end run. (All of the committed test data together comes to under 500 KB, so it is cheap to keep in the repository.)

```bash
# Run from the repo root (inside Singularity container or with all tools on PATH)
./tests/run_tests.sh all         # all layers
./tests/run_tests.sh schema      # schema validation only (fastest; no biology tools needed)
./tests/run_tests.sh unit        # R-script unit tests
```

See [`tests/README.md`](tests/README.md) for the full instructions, including how to regenerate the fixtures.

## Required config entries

The entries below are the ones worth understanding before your first run.

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

### `filtering / tss_srna_sizes` and the read-size composition filter

Abundant uncapped small RNAs are the one contaminant the other filters cannot
reach. In plants the 21-25 nt siRNAs, above all the 24 nt Pol IV class over
transposons, survive TEX/AP well enough to be called as TSS clusters, and
enrichment over the input cannot reject them wherever the input library has too
little coverage to measure a local background. Their read length gives them
away, because genuine initiation is not confined to one small size class.

List the contaminating read lengths in `tss_srna_sizes`, and a TSS is dropped
from `tss.final.bed` when more than `tss_max_srna_fraction` of its reads fall in
those sizes. Only libraries with at least `tss_srna_min_reads` reads in the
cluster get a vote, and `tss_srna_min_samples` of them must agree before the
cluster goes. The default `[]` disables the filter, and while it is disabled the
`tss_size_composition` rule does not read the tag directories at all.

```yaml
filtering:
  tss_srna_sizes: [21, 22, 23, 24, 25]
  tss_max_srna_fraction: 0.5
  tss_srna_min_reads: 100
  tss_srna_min_samples: 1
```

The default `tss_srna_min_samples: 1` drops a cluster that any single library
flags, which is the sensitive choice when libraries differ in size selection: a
40-80 nt library cannot see a 24 nt population that a 20-80 nt library shows
plainly, so demanding agreement would let the contaminant through. Raise it if
you would rather keep a cluster that only one library objects to.

Per-cluster counts land in `tss.consensus.sizes.txt` whether or not anything is
removed, so the size composition can be inspected directly.

### `filtering / tss_max_top_sizes_fraction` and the top-lengths filter

The size list above only catches contaminants whose length you can name in advance. This
second filter needs no such list: it asks how concentrated a cluster is in its own
commonest read lengths. Genuine initiation is heterogeneous, because a promoter fires
across a window and the resulting RNAs vary in length, so a protein-coding cluster spreads
over tens of lengths and its two commonest hold only about a fifth of it. A discretely
processed RNA has one 5′ end and one 3′ end and puts nearly everything into one or two
lengths.

Set `tss_max_top_sizes_fraction` below 1 to enable, and `tss_top_sizes_n` for how many
lengths to add up. The same `tss_srna_min_reads` and `tss_srna_min_samples` gates apply.

```yaml
filtering:
  tss_top_sizes_n: 2
  tss_max_top_sizes_fraction: 0.8
```

Measured on Arabidopsis csRNA-seq (113,557 consensus clusters, libraries with at least 30
reads in a cluster), those settings flag **0.09% of protein-coding clusters and 17% of
transposon clusters**. The medians it works from: pcTSS clusters spread over 33 distinct
lengths with a top-2 share of 0.21, while transposon clusters manage 10 lengths and 0.49.
It catches, for example, an abundant 27-28 nt species over a BRODYAGA1A element whose 5′
and 3′ ends are both fixed (top-2 share 0.83), which sits outside any plant small RNA size
class and so escapes `tss_srna_sizes: [21, 22, 23, 24, 25]` entirely.

Concentration is estimated from few reads at weakly expressed clusters, so the flag rate is
depth-dependent: pooled across classes it runs 2.4% at 30-50 reads per cluster, 1.1% at
51-100, and 0.1-0.3% above 300. That is why `tss_srna_min_reads` defaults to 100. Real
single-locus contaminants are usually very abundant (the BRODYAGA1A species carries over
20,000 reads in one library), so the floor keeps them while dropping most of the
small-sample noise. Lower it to around 30 if you would rather be sensitive.

Note that `tss_top_sizes_n` is baked into `tss.consensus.sizes.txt` when it is written, so
changing it re-runs the `tss_size_composition` rule.

## Workflow Steps

| Rule | Description |
|------|-------------|
| `build_index` | (Conditional) Build a genome index for the selected aligner from `genome_fasta` if no index is found at `genome_index`. |
| `trim` | Adapter trimming and length truncation (bfqutils). Stats written to `qc/`. |
| `align` | Genome alignment (bwa-aln / bwa-mem / STAR / bowtie2 / hisat2), MAPQ and length filtering, BAM sorting and indexing. |
| `make_tagdir` | Tag directory from BAM (bam2td). |
| `find_tss_initial` | Per-sample TSS calling (HOMER `findcsRNATSS.pl`) using the paired input library as background. Run on both csRNA and input samples. |
| `make_raw_bedgraph` | Strand-separated raw bedGraph files (HOMER `makeUCSCfile`). |
| `gather_stats` | Extract total reads, organellar reads, and tag frequencies from HOMER tagInfo files. |
| `merge_initial_tss` | Strand-aware merge of all per-sample TSS BEDs. Chromosomes not in `chrom_sizes` are excluded. |
| `quantify_initial_{cs,in}_tss` | HOMER `annotatePeaks.pl` quantification of merged TSS sets in all libraries. |
| `qc_initial_tss` | Calculate FRiP, nuclear read %, csRNA enrichment, miRNA depletion, and phosphorylation efficiency. Samples with FRiP < `min_cs_frip` or nuclear % < `min_pct_nuclear` are flagged as FAIL. |
| `collect_consensus_tss` | Build the consensus TSS set (`tss.consensus.bed`) from csRNA samples, requiring detection in at least `tss_min_reps` replicates. TSSs narrower than 150 bp are padded; overlapping TSSs are split. TSSs overlapping miRNA / pre-tRNA loci (if `qc / mirnas` and `qc / trnas` are set) are removed. |
| `quantify_final_tss` | HOMER quantification of the consensus TSSs in all csRNA libraries. |
| `tss_size_composition` | Per-cluster read counts, small-RNA-sized read counts and top-lengths concentration from the tag directories, written to `tss.consensus.sizes.txt`. Reads nothing when both size filters are off. |
| `normalize_tss_quantification` | CPM filter (drop TSSs below `filtering / tss_min_cpm` in fewer than `filtering / tss_min_samples` csRNA samples), then the two read-size composition filters (`filtering / tss_srna_sizes` and `filtering / tss_max_top_sizes_fraction`), then TMM normalization (edgeR `TMMwsp`) of the retained counts. Writes the filtered set as `tss.final.bed`. All of these default to no-ops, keeping every TSS. |
| `quantify_final_in_tss` | HOMER quantification of the filtered `tss.final.bed` in all input libraries. Used by the final QC step. |
| `qc_final_tss` | Compute FRiP, csEnrichment, replicate correlations, and TSS detection rate on the filtered `tss.final.bed`. Writes `qc/qc_final_cs.txt` and `qc/qc_final_in.txt`. |
| `generate_normalized_bw` | Multiply raw bedGraphs by RPM scale factors and export as bigWig. Small RNA regions can optionally be masked. |

## Output Files

All outputs land in `files / output_dir` (default: `results/`):

| File | Description |
|------|-------------|
| `tss.consensus.bed` | Unfiltered consensus TSS set (before CPM filter). |
| `tss.consensus.sizes.txt` | Per csRNA library and cluster: total reads (`.reads`), reads in `tss_srna_sizes` (`.srna`), and reads in the cluster's `tss_top_sizes_n` commonest lengths (`.topn`). Cluster names only when both size filters are off. |
| `tss.final.bed` | Filtered consensus TSS coordinates (BED6). With default `tss_min_cpm: 0` and `tss_srna_sizes: []` this equals `tss.consensus.bed`. |
| `tss.final.raw.txt` | Raw tag counts per TSS per csRNA sample (filtered set only). |
| `tss.final.cpm.txt` | TMM-normalized CPM counts. |
| `norm_factors.txt` | edgeR TMM normalization factors and RPM multipliers. |
| `bw/{sample}.rpm.pos.bw` | Forward-strand RPM-normalized bigWig. |
| `bw/{sample}.rpm.neg.bw` | Reverse-strand RPM-normalized bigWig (scores are negative). |
| `qc/qc_initial_cs.txt` | Per-sample QC metrics for csRNA libraries computed on the merged initial TSS set (FRiP, enrichment, etc.). |
| `qc/qc_initial_in.txt` | Per-sample QC metrics for input libraries computed on the merged initial TSS set. |
| `qc/qc_final_cs.txt` | Per-sample QC metrics for csRNA libraries re-derived from the filtered `tss.final.bed`. Contains the same columns as `qc_initial_cs.txt` plus the `Final*` and `NConsensusTSS`/`NFinalTSS`/`NFilteredTSS` columns. All metrics (FRiP, enrichment, contamination) are re-computed from raw counts, not copied from the initial QC. |
| `qc/qc_final_in.txt` | Per-sample QC metrics for input libraries, augmented with `Final*` columns. |
| `qc/replicate_correlation.txt` | Pairwise Spearman and Pearson correlations between csRNA replicates of the same `sample_name`, for both the initial TSS set (Stage=Initial) and the filtered set (Stage=Final). Empty if no `sample_name` has ≥ 2 replicates. |
| `qc/stats_initial_cs.txt` | Raw alignment and tag-count statistics for csRNA libraries (TotalReads, OrganelleReads, PosReads, NegReads). |
| `qc/stats_initial_in.txt` | Raw alignment and tag-count statistics for input libraries. |
| `qc/{sample}.trimming.txt` | bfqutils trimming summary. |
| `qc/{sample}.aln.txt` | samtools flagstat alignment summary. |
| `qc/{sample}.aligner.log` | Aligner stderr (bwa/STAR/bowtie2/hisat2). |

## QC Metrics

The `qc_initial_cs.txt` / `qc_initial_in.txt` tables contain the columns below. `qc_final_cs.txt` / `qc_final_in.txt` contain the same columns plus the `Final*` / `NConsensusTSS` / `NFinalTSS` / `NFilteredTSS` columns, with all values re-derived from the filtered TSS set rather than copied from the initial QC:

| Column | Description |
|--------|-------------|
| `csFRiP` | Fraction of nuclear reads falling in csRNA-called TSSs. Should be > 0.9 for a good csRNA library. |
| `sFRiP` | Fraction of nuclear reads in input-called (small RNA) TSSs. |
| `PctNuclear` | Percentage of reads mapping to nuclear chromosomes (i.e. excluding `qc / organelle_chroms`). |
| `csEnrichment` | Ratio of csRNA FRiP to input FRiP, i.e. how enriched the capped initiation signal is relative to the background. |
| `sDepletion` | Inverse ratio of small RNA signal between csRNA and input libraries. |
| `PretRNAPct` | Fraction of input-TSS reads overlapping pre-tRNA loci (requires `qc / trnas`). |
| `PhosEfficiency` | Ratio of pre-tRNA depletion in csRNA vs input; a measure of 5′-phosphate removal efficiency. |
| `miRNADepletion` | Ratio of miRNA depletion in csRNA vs input, an independent phosphorylation-efficiency metric (requires `qc / mirnas`). |
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
