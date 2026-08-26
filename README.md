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

### `program / homer / tss / gtf` (optional, and worth setting)

Optional gene annotation passed to HOMER as `-gtf`. Without it, and without a HOMER genome
that carries its own annotation, HOMER cannot build the annotated-TSS and exon sets it uses
as true and false positives to choose the enrichment threshold for each library. It falls
back to `default_log2_fold` silently: the only sign in the log is

```
Total TP (tss) regions: 0     Total FP (exon) regions: 0
Maximum CDF difference: -10000000000
Using the default input threshold (1)
```

and, in the per-library `tss/<sample>.stats.txt`, `Fraction Promoter-Distal TSS
clusters: 100.00%` with `Fraction of stable transcript TSS clusters: 0.00%`, which are
placeholders rather than results. The pipeline now rewrites both to `na` when HOMER had
nothing to compute them from, so the file no longer reports a placeholder as a
measurement; see the note under Output Files.

With a GTF, the same library reports non-empty sets and a threshold chosen from the data.
On the synthetic test genome:

```
Custom annotation GTF file: .../annotation.gtf (using transcript_id)
Total TP (tss) regions: 15     Total FP (exon) regions: 5
Maximum CDF difference: 0.8
log2 fold vs. input: 0.415031488137169
```

Ignore the `Skipping TSS assignment (can't find file for genome ...)` line while judging
this: it comes from the `annotatePeaks` calls inside HOMER wanting a full HOMER genome
directory for their own TSS-distance annotation, and it appears whether or not `-gtf` was
given. The TP and FP counts are the signal.

Supplying a GTF is the largest single quality lever in the configuration, and it **changes
your TSS calls**, so compare a GTF run against your existing consensus before adopting it
rather than swapping mid-project.

```yaml
program:
  homer:
    tss:
      gtf: "/path/to/annotation.gtf"
```

Empty by default, which is the pipeline's previous behaviour.

### `program / homer / tss / program`, `pseudo_count`, `default_log2_fold`, `local_fold`

`findcsRNATSS.pl` prints *"this program is a legacy placeholder - please use
findcsRNATSR.pl instead"*. The two are the same program: 81 lines of diff, all renames and
help text. Verified on a 34,891-cluster library, their tables are identical row for row, so
switching changes nothing except output file names and labels, which the pipeline
normalises. The default stays on the legacy name because `findcsRNATSR.pl` ships with an
unpatched `use lib` pointing at its author's own install, so it only finds `HomerConfig`
when HOMER's `bin` is on `PERL5LIB`. Check that inside your container, then set
`program: "findcsRNATSR.pl"`.

Three HOMER thresholds are now configurable, all defaulting to HOMER's own values so
nothing changes until you touch them:

- `pseudo_count` (1.0) is added to both sides of the enrichment ratio, and is what lets a
  cluster with zero input coverage pass automatically, since `log2((cs+1)/(0+1))` is large
  by construction. Raise it to demand real input coverage.
- `default_log2_fold` (1) is the threshold actually in force whenever no annotation is
  supplied, per the section above.
- `local_fold` (2) is the local enrichment required during initial TSS identification.

### `filtering / exclude_failed_from_consensus`

The consensus set is a union, so a library that failed QC does not simply have poor numbers
of its own: the spurious clusters it calls enter the shared set, and everything downstream
quantifies them. On one Arabidopsis panel, two contaminated libraries contributed about
2,300 transposon-promoter clusters each, and 1,745 of 3,839 such clusters came from exactly
one library.

Default `False`, which keeps the historical behaviour of building the union from every
csRNA library. A warning naming the offending libraries is printed either way.

```yaml
filtering:
  exclude_failed_from_consensus: true
```

`Status` is FAIL when `csFRiP` is below `qc / min_cs_frip`, nuclear percentage is below
`qc / min_pct_nuclear`, or the enrichment threshold HOMER chose is below
`qc / min_log2_fold`; `StatusReason` says which. Consider `tss_min_reps: 2` alongside it,
so that no single library can inject a cluster on its own.

### `qc / min_log2_fold`

Minimum enrichment threshold a library's TSS calling may use before the library fails QC.
Unset by default, in which case it falls back to `program / homer / tss / default_log2_fold`.

With an annotation, HOMER chooses this threshold per library from the data rather than
using `default_log2_fold`, and writes its choice into `tss/<sample>.stats.txt` as
`log2 fold vs. input`. For sound libraries it picks sensible values (log2 1.4 to 1.9 on one
Arabidopsis panel) and cluster counts drop 4 to 10%. For libraries that have lost their
capped signal it picks thresholds at or below zero:

```
library                chosen log2FC    valid clusters, no GTF -> with GTF
polyphos_tex_csrna1            0.144     26,025 -> 43,454   (+67%)
polyphos_tex_csrna2            0.007      7,613 -> 28,816  (+278%)
tex_only_csrna1               -0.680     14,229 -> 28,864  (+103%)
tex_only_csrna2               -0.056     18,874 -> 24,677   (+31%)
```

A negative threshold accepts clusters carrying less signal than their own input, which
cannot be a TSS under any reading, and the annotation therefore makes a bad library look
more productive rather than less. Because the consensus is a union, those calls propagate:
on that panel 25,871 of 80,887 final clusters rested on those four libraries alone. The
threshold is chosen inside HOMER and there is no flag to floor it, so the pipeline reads it
back and gates on it.

The fallback is the defensible bound rather than zero: a data-driven threshold should never
end up looser than the constant it replaced. Set the key explicitly to gate independently
of `default_log2_fold`.

```yaml
qc:
  min_log2_fold: 1.0
```

A library that trips this gate is FAIL with `Log2FoldThreshold` in its `StatusReason`, and
`filtering / exclude_failed_from_consensus` is what keeps its clusters out of the union.

### `program / keep_unmapped_sample`

Number of unmapped reads to keep per library, as `qc/{sample}.unmapped.sample.fastq.gz`.
Default 0 writes nothing.

Two thirds of a csRNA-seq library can be lost at alignment (a median 74.9% per library on
one Arabidopsis panel), and with `keep_bam`, `keep_trimmed_fastq` and
`bwa_aln / keep_sai` all false, neither the aligned nor the unaligned reads survive the
run, so a contaminant screen or a low-complexity check afterwards means aligning again.
A few hundred thousand reads costs a few MB and answers the question.

```yaml
program:
  keep_unmapped_sample: 200000
```

`qc/{sample}.aln.raw.txt` says how large that loss is and how it splits between reads that
never aligned, reads below the MAPQ floor and reads dropped by the flag and length
filters; the same three numbers appear in the QC tables as `AlignedReads`,
`BelowMapqReads` and `FilteredOutReads`. This sample is for asking what the unaligned
reads actually are.

### `qc / min_cs_frip_final` and `qc / min_pct_nuclear_final`

The initial and final QC tables both report `csFRiP`, but they measure it against
different regions: the initial merged TSS set in one, `tss.final.bed` in the other. One
number therefore cannot gate both. Every library's `csFRiP` falls when the consensus set
shrinks (0.030 to 0.058 on one 22-library Arabidopsis panel, purely from the set getting
smaller), so a gate tuned so that the initial table flags exactly the four collapsed
libraries went on to flag 8 of 22 in the final table, four of them sound.

Both keys are unset by default, in which case they fall back to `min_cs_frip` and
`min_pct_nuclear`, so an existing config behaves exactly as it did before they existed.

```yaml
qc:
  min_cs_frip: 0.80         # gates the initial table, and the consensus through it
  min_cs_frip_final: 0.70   # gates the final table only
```

Only the initial `Status` feeds `collect_consensus_tss`; see the note at the end of the QC
section.

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
transposons, survive the enzymatic depletion well enough to be called as TSS
clusters. Enrichment over the input holds them off while the csRNA library is
clean, but it degrades as the library degrades, because a library that has lost
its capped signal is proportionally richer in siRNA than its own input is.
Measured on Arabidopsis, the share of siRNA-dominated clusters passing a 2-fold
enrichment filter runs 1% in a good library, 7% without the AP phosphatase, and
61% in the worst library in the panel. Read length is the signal that does not
degrade with the library, because genuine initiation is never confined to one
small size class.

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
| `align` | Genome alignment (bwa-aln / bwa-mem / STAR / bowtie2 / hisat2), MAPQ and length filtering, BAM sorting and indexing. Counts the unfiltered alignment on the way past into `qc/{sample}.aln.raw.txt`, which is where the reads lost at this step are accounted for. |
| `make_tagdir` | Tag directory from BAM (bam2td). |
| `find_tss_initial` | Per-sample TSS calling (HOMER `findcsRNATSS.pl`) using the paired input library as background. Run on both csRNA and input samples. |
| `make_raw_bedgraph` | Strand-separated raw bedGraph files (HOMER `makeUCSCfile`). |
| `gather_stats` | Extract total reads, organellar reads, and tag frequencies from HOMER tagInfo files, plus the alignment loss (`AlignedReads`, `BelowMapqReads`, `FilteredOutReads`) from the per-library alignment summaries. |
| `merge_initial_tss` | Strand-aware merge of all per-sample TSS BEDs. Chromosomes not in `chrom_sizes` are excluded. |
| `quantify_initial_{cs,in}_tss` | HOMER `annotatePeaks.pl` quantification of merged TSS sets in all libraries. |
| `qc_initial_tss` | Calculate FRiP, nuclear read %, csRNA enrichment, miRNA depletion, and phosphorylation efficiency. Samples with FRiP < `min_cs_frip`, nuclear % < `min_pct_nuclear`, or an enrichment threshold below `min_log2_fold` are flagged as FAIL, with `StatusReason` naming the gates that tripped. The gates reach the script as rule params, so editing one and rerunning with `--rerun-triggers params` (Snakemake's default set includes `params`) recomputes the table. |
| `collect_consensus_tss` | Build the consensus TSS set (`tss.consensus.bed`) from csRNA samples, requiring detection in at least `tss_min_reps` replicates. TSSs narrower than 150 bp are padded; overlapping TSSs are split. TSSs overlapping miRNA / pre-tRNA loci (if `qc / mirnas` and `qc / trnas` are set) are removed. |
| `quantify_final_tss` | HOMER quantification of the consensus TSSs in all csRNA libraries. |
| `tss_size_composition_library` | Per-cluster read counts, small-RNA-sized read counts and top-1..N length concentration for one csRNA library, from its tag directory. One job per library, so Snakemake runs them in parallel. Reads nothing when both size filters are off. |
| `tss_size_composition` | Joins the per-library tables into `tss.consensus.sizes.txt`. |
| `normalize_tss_quantification` | CPM filter (drop TSSs below `filtering / tss_min_cpm` in fewer than `filtering / tss_min_samples` csRNA samples), then the two read-size composition filters (`filtering / tss_srna_sizes` and `filtering / tss_max_top_sizes_fraction`), then TMM normalization (edgeR `TMMwsp`) of the retained counts. Writes the filtered set as `tss.final.bed`. All of these default to no-ops, keeping every TSS. |
| `quantify_final_in_tss` | HOMER quantification of the filtered `tss.final.bed` in all input libraries. Used by the final QC step. |
| `qc_final_tss` | Compute FRiP, csEnrichment, replicate correlations, and TSS detection rate on the filtered `tss.final.bed`. Gated by `qc / min_cs_frip_final` and `qc / min_pct_nuclear_final`, each falling back to its initial counterpart. Writes `qc_final_cs.txt` and `qc_final_in.txt`. |
| `generate_normalized_bw` | Multiply raw bedGraphs by RPM scale factors and export as bigWig. Small RNA regions can optionally be masked. |
| `run_info` | Write `run_info.txt`: pipeline version and git commit, host, resolved config, sample-table checksum, and the version of every tool the run used. |

## Output Files

All outputs land in `files / output_dir` (default: `results/`):

| File | Description |
|------|-------------|
| `tss.consensus.bed` | Unfiltered consensus TSS set (before CPM filter). |
| `tss.consensus.sizes.txt` | Per csRNA library and cluster: total reads (`.reads`), reads in `tss_srna_sizes` (`.srna`), and reads in the cluster's 1..`tss_top_sizes_max` commonest lengths (`.top1` … `.topN`). All n are precomputed so `tss_top_sizes_n` can be retuned without re-reading the tag directories. Cluster names only when both size filters are off. |
| `run_info.txt` | Provenance for the run: pipeline version and commit, host, tool versions, sample-table checksum and the fully resolved config. |
| `tss.final.bed` | Filtered consensus TSS coordinates (BED6). With default `tss_min_cpm: 0` and `tss_srna_sizes: []` this equals `tss.consensus.bed`. |
| `tss.final.raw.txt` | Raw tag counts per TSS per csRNA sample (filtered set only). |
| `tss.final.cpm.txt` | TMM-normalized CPM counts. |
| `norm_factors.txt` | edgeR TMM normalization factors and RPM multipliers. |
| `bw/{sample}.rpm.pos.bw` | Forward-strand RPM-normalized bigWig. |
| `bw/{sample}.rpm.neg.bw` | Reverse-strand RPM-normalized bigWig (scores are negative). |
| `qc/qc_initial_cs.txt` | Per-sample QC metrics for csRNA libraries computed on the merged initial TSS set (FRiP, enrichment, etc.). |
| `qc/qc_initial_in.txt` | Per-sample QC metrics for input libraries computed on the merged initial TSS set. |
| `qc_final_cs.txt` | Per-sample QC metrics for csRNA libraries re-derived from the filtered `tss.final.bed`. Same columns as `qc_initial_cs.txt` plus `FinalTSSDetected`, `NConsensusTSS`, `NFinalTSS` and `NFilteredTSS`. All metrics are re-computed from raw counts, not copied from the initial QC, so `csRiP` / `csFRiP` / `csEnrichment` here refer to the **final** TSS set while the same names in `qc_initial_cs.txt` refer to the initial merged set. Do not compare the two files column-by-column. |
| `qc_final_in.txt` | Per-sample QC metrics for input libraries, augmented with `Final*` columns. |
| `qc/replicate_correlation.txt` | Pairwise Spearman and Pearson correlations between csRNA replicates of the same `sample_name`, for both the initial TSS set (Stage=Initial) and the filtered set (Stage=Final). Empty if no `sample_name` has ≥ 2 replicates. |
| `qc/stats_initial_cs.txt` | Raw alignment and tag-count statistics for csRNA libraries (RawReads, TrimmedReads, AlignedReads, BelowMapqReads, FilteredOutReads, TotalReads, OrganelleReads, PosReads, NegReads). |
| `qc/stats_initial_in.txt` | Raw alignment and tag-count statistics for input libraries. |
| `tss/{sample}.stats.txt` | HOMER's own per-library report from TSS calling: cluster counts, the enrichment threshold it chose (`log2 fold vs. input`), and the true/false-positive set sizes it chose it from. Post-processed by the pipeline: the promoter-distal and stable-transcript fractions are rewritten to `na` when HOMER had no annotation and no RNA-seq to compute them from, since it otherwise reports its placeholders (100.00% and 0.00%) as though they were measurements. Anything parsing this file should expect `na` in those two positions. |
| `qc/{sample}.trimming.txt` | bfqutils trimming summary. |
| `qc/{sample}.aln.txt` | samtools flagstat, run on the **filtered** BAM. It therefore reports 100% mapped for every library and says nothing about what the filter removed; `aln.raw.txt` is the file that does. |
| `qc/{sample}.aln.raw.txt` | Accounting of the unfiltered alignment: reads seen, unmapped, secondary, supplementary, reads with a primary alignment, then how many of those fell below the MAPQ floor, failed the remaining filters, and passed. Ends with a MAPQ histogram and the filter settings in force. Counted in the pipe, since the unfiltered alignment is never written to disk. |
| `qc/{sample}.unmapped.sample.fastq.gz` | A sample of unmapped reads, only when `program / keep_unmapped_sample` is above zero. |
| `qc/{sample}.aligner.log` | Aligner stderr (bwa/STAR/bowtie2/hisat2). For STAR, `Log.final.out` is appended to it. |

## QC Metrics

The `qc_initial_cs.txt` / `qc_initial_in.txt` tables contain the columns below. `qc_final_cs.txt` / `qc_final_in.txt` contain the same columns plus `FinalTSSDetected` / `NConsensusTSS` / `NFinalTSS` / `NFilteredTSS`, with all values re-derived from the filtered TSS set rather than copied from the initial QC. That means a column such as `csFRiP` measures the initial merged TSS set in the initial tables and the final set in the final tables, so the two files are not comparable column-by-column:

| Column | Description |
|--------|-------------|
| `Log2FoldThreshold` | The enrichment threshold this library's TSS calling actually used, read back from `tss/<sample>.stats.txt`. Chosen by HOMER from the data when an annotation is available, otherwise `program / homer / tss / default_log2_fold`. Gated by `qc / min_log2_fold`. `NA` when HOMER wrote no threshold at all, which is what happens for a library with no valid clusters: it hits a division by zero on the way to the promoter-distal fraction and its report stops early. `NA` counts as a failure, and input libraries commonly show it. |
| `AlignedReads` | Reads with a primary alignment, counted before filtering. `TrimmedReads` minus this is the reads that never aligned. |
| `BelowMapqReads` | Aligned reads discarded by `filtering / alignment_mapq`, which in practice means multi-mappers. |
| `FilteredOutReads` | Aligned reads that cleared MAPQ and were then discarded by `exclude_flags`, `include_flags`, `min_alignment_length` or `max_mismatch`. |
| `csFRiP` | Fraction of nuclear reads falling in csRNA-called TSSs. Should be > 0.9 for a good csRNA library. |
| `sFRiP` | Fraction of nuclear reads in input-called (small RNA) TSSs. |
| `PctNuclear` | Percentage of reads mapping to nuclear chromosomes (i.e. excluding `qc / organelle_chroms`). |
| `csEnrichment` | Ratio of csRNA FRiP to input FRiP, i.e. how enriched the capped initiation signal is relative to the background. |
| `sDepletion` | Inverse ratio of small RNA signal between csRNA and input libraries. |
| `PretRNAPct` | Pre-tRNA reads as a percentage of pre-tRNA reads plus reads in TSSs, i.e. contaminant per unit of signal rather than per library read (requires `qc / trnas`). |
| `PhosEfficiency` | Ratio of pre-tRNA depletion in csRNA vs input; a measure of 5′-phosphate removal efficiency. |
| `miRNADepletion` | Ratio of miRNA depletion in csRNA vs input, an independent phosphorylation-efficiency metric (requires `qc / mirnas`). |
| `miRNAPct` | miRNA reads as a percentage of miRNA reads plus reads in TSSs. Same basis as `PretRNAPct`, so the two contamination percentages are directly comparable. |
| `StrandBalance` | Fraction of mapped reads on the + strand. Expect ~0.5 in csRNA libraries; large deviations flag adapter contamination, library-prep strand bias, or pile-ups at a few highly expressed loci. Looser bounds in input libraries since small-RNA biology is genuinely strand-skewed. |
| `TSSDetected` | Fraction of merged-set TSSs with ≥ 1 tag in this library. Low values flag undersequenced libraries. |
| `FinalTSSDetected` | (`qc_final_*` only) Fraction of `tss.final.bed` TSSs with ≥ 1 tag in this library. |
| `NConsensusTSS` | (`qc_final_*` only) Number of TSSs in `tss.consensus.bed` (before CPM filter). |
| `NFinalTSS` | (`qc_final_*` only) Number of TSSs retained in `tss.final.bed` after CPM filter. |
| `NFilteredTSS` | (`qc_final_*` only) Number of TSSs dropped by the CPM filter. |
| `MinReplCorrSpearman` | Minimum Spearman correlation between this sample's TSS counts and any other csRNA replicate of the same `sample_name`. `NA` when only one replicate exists. Sharp drops (e.g. < 0.9) flag sample swaps or replicate dropouts. csRNA samples only. |
| `MinReplCorrPearson` | Same as above but Pearson correlation on log1p-transformed counts. |
| `Status` | `Ok` or `FAIL`, from the gates listed under `StatusReason`. A gate that cannot be evaluated (an `NA` csFRiP, say) counts as a failure rather than as a pass. csRNA samples only. |
| `StatusReason` | Comma-separated list of the gates that tripped, empty when `Status` is `Ok`. This is what makes a FAIL actionable without re-deriving every gate by hand. |
| `StatusBasis` | `initial` in `qc_initial_cs.txt` and `final-advisory` in `qc_final_cs.txt`. `collect_consensus_tss` reads the **initial** `Status` and nothing reads the final one, so a FAIL in the final table has already had no effect on the TSS set that table describes. |

### Reading these numbers, and what they do not mean

Four things about this table mislead people, including the person who wrote most of it.

**The same column name means different regions in the two files.** `csRiP`, `csFRiP` and
`csEnrichment` are computed against the initial merged TSS set in `qc_initial_*.txt` and
against `tss.final.bed` in `qc_final_*.txt`. Both are correct; they are just not the same
measurement, so never join or compare the two files column-by-column.

**The contamination percentages are signal-relative, not library-relative.** `PretRNAPct`
and `miRNAPct` are both contaminant reads over contaminant plus reads-in-TSS, so they answer
"per read I care about, how much junk", not "what fraction of the library is junk". The
library-fraction version is smaller by roughly the csFRiP factor, which is itself lower in
worse libraries, so a contaminated library looks somewhat worse on these than a
library-fraction metric would show.

**A high duplicate rate is not automatically bad, and a low one is not automatically good.**
Duplication at matched depth is complexity convolved with concentration: reads pile up where
the signal is, so a library that puts more of its reads into TSSs duplicates more at the same
depth. On one Arabidopsis panel the *least* duplicated library was among the two worst by
every other measure, because its reads were scattered over degradation background. If you
want complexity, compare distinct positions at a matched depth restricted to real promoters.

**Every per-class figure counts uniquely-mapping reads only.** All of them are computed
from the filtered BAM, which is past the MAPQ floor, so `PretRNAPct`, `miRNAPct` and any
rRNA figure measure the unique-mapping tail of those species rather than their abundance.
On one Arabidopsis panel rRNA came out at 0.07% of aligned 5′ ends, which is not credible
for a total-RNA-derived library and is better read as evidence that the rRNA is in the
discarded multi-mapping fraction: TAIR12 assembles the rDNA arrays as roughly 4,800 rRNA
gene entries, which no short read can place uniquely. `BelowMapqReads` is how much of the
library that fraction is.

**`csEnrichment` is not a reliable discriminator at these group sizes.** On a panel of two to
three libraries per condition it could not separate 1.97 from 2.00 for a treatment that
demonstrably changed the libraries. The metrics that did separate consistently were
`PretRNAPct`, `PhosEfficiency`, alignment survival and the useful-read fraction. Power
comparisons on those.

**Only the initial `Status` acts on anything.** `collect_consensus_tss` reads
`qc_initial_cs.txt`, so that is the table whose gates (`qc / min_cs_frip`,
`qc / min_pct_nuclear`) can change a run, and `filtering / exclude_failed_from_consensus`
is what makes them do so. The final table's `Status` is advice, labelled as such in
`StatusBasis`, and it has its own gates because it measures a different thing.
