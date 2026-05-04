# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- End-to-end tests for bwa-mem and hisat2 alignment programs, runnable via
  `./tests/run_tests.sh bwa-mem` and `./tests/run_tests.sh hisat2`; both
  are included in the `all` subcommand.

## [0.1.0] - 2026-05-04

### Added
- Snakemake workflow with 13 rules covering the full csRNA-seq analysis:
  adapter trimming (bfqutils), genome alignment (bwa-aln, bwa-mem, STAR,
  bowtie2, hisat2), HOMER tag directory creation, initial TSS detection,
  strand-separated bedGraph generation, alignment statistics, TSS merging,
  quantification, TMM normalisation, consensus TSS collection, CPM filtering,
  post-filter QC, and normalised bigWig generation.
- QC metrics: FRiP, nuclear read fraction, csRNA enrichment over input,
  miRNA/tRNA contamination, Pearson replicate correlation, and per-sample
  alignment statistics.
- `input_name` column in the sample table to allow multiple conditions to
  share a single input library.
- CPM filter (`tss_min_cpm`, `tss_min_samples`) to remove low-confidence TSSs
  from the consensus set, with a dedicated post-filter QC step.
- Small-RNA masking option (`mask_srnas_from_bw`) to exclude miRNA/tRNA loci
  from normalised bigWig tracks.
- Schema validation for the sample table and configuration file via
  `snakemake.utils.validate`.
- `conda` environment definition (`workflow/envs/env.yaml`).
- Test suite under `tests/` with four layers: schema validation, Snakemake
  dry-run validation, R-script unit tests, and end-to-end tests (bwa-aln,
  STAR, bowtie2) against a synthetic 80 kb genome.

### Fixed
- STAR alignment: bfqutils appends a trailing space to every read ID line;
  STAR 2.7.x rejects such files.  Reads are now piped through `sed` inside a
  bash process substitution to strip trailing whitespace before STAR sees them.
- STAR output prefix derivation: using `Path.with_suffix()` on a Snakemake
  placeholder string corrupted the closing `}`.  The prefix is now computed at
  shell runtime via `dirname`/`basename`.
- `generate_normalized_bw.R`: negative-strand bedGraphs are no longer empty
  when all synthetic reads are forward-strand, preventing a crash in
  `rtracklayer::export.bw()`.
- `zcat` on macOS (BSD) cannot decompress `.gz` files; replaced with
  `gzip -dc` throughout.

[0.1.0]: https://github.com/bjmt/csRNAseq_pipeline/releases/tag/v0.1.0
