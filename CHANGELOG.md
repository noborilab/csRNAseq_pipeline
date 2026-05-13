# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.4.0] - 2026-05-13

### Changed
- `make_tagdir` now uses `bam2td` instead of HOMER `makeTagDirectory`.
  bam2td reads sorted BAMs directly (no intermediate SAM conversion) and
  produces a compatible HOMER tag directory.  The `shadow: "minimal"` and
  the `samtools view -h` conversion step are removed.  `makeTagDirectory -update -checkGC -genome` is run afterwards to
  populate QC files (e.g. `tagFreq.txt`, `tagGCcontent.txt`) that bam2td
  does not generate.  The `homer.tagdir.extra_args` config key is still
  supported and passed directly to bam2td.

## [0.3.0] - 2026-05-13

### Added
- `RawReads` and `TrimmedReads` columns in `stats_initial_cs.txt` and
  `stats_initial_in.txt` (and propagated to all QC output files), parsed
  from bfqutils `.trimming.txt` logs; multi-file samples are handled by
  summing across all per-file log entries.
- `miRNAPct` column (miRNA reads / NuclearReads × 100) in all QC output
  files, computed alongside the existing `miRNA` and `miRNADepletion` metrics.
- `tss_min_cs_in_ratio` and `tss_min_ratio_samples` config parameters: a
  new csRNA/input read-ratio filter applied in `normalize_tss_quantification`
  after the CPM filter.  Uses the already-computed input quantification of
  the initial csRNA TSS set (`all_cs.tss_merged_quant_in.txt`) so no extra
  HOMER run is required.  Defaults to 0 (disabled).

### Changed
- `key` column renamed to `group` in all QC output files (`qc_initial_cs.txt`,
  `qc_initial_in.txt`, `qc_final_cs.txt`, `qc_final_in.txt`).
- `qc_final_cs.txt` and `qc_final_in.txt` now written to the results root
  directory rather than `qc/`.
- `tss_merged/` outputs (`all_cs.tss_merged.bed`, `all_in.tss_merged.bed`,
  and all derived quantification files) now land in the results directory
  rather than the intermediate directory.
- `generate_normalized_bw` rule memory increased from 2 000 MB to 4 000 MB.

## [0.2.0] - 2026-05-06

### Added
- End-to-end tests for bwa-mem and hisat2 alignment programs, runnable via
  `./tests/run_tests.sh bwa-mem` and `./tests/run_tests.sh hisat2`; both
  are included in the `all` subcommand.
- `workflow/envs/minimal_env.yaml` rewritten as a portable, cross-platform
  conda spec with correct channel order (`conda-forge` before `bioconda`),
  minimum-version constraints, and notes for the two manual-install
  dependencies (bfqutils, HOMER on macOS).
- `workflow/envs/full_env.yaml` for development and running the full test
  suite; adds STAR, bowtie2, hisat2, pytest, snakefmt, and vim on top of
  the minimal spec.
- `qc_final_cs.txt` and `qc_final_in.txt` now contain the full set of QC
  metrics (csRiP, sRiP, csFRiP, sFRiP, csRNACappedPct, csEnrichment,
  sDepletion, PretRNA, PretRNAPct, PhosEfficiency, miRNA, miRNADepletion)
  re-derived from raw quantification files rather than copied from the
  initial QC step, so they reflect the CPM-filtered TSS set.
- `qc/replicate_correlation.txt` now contains pairwise correlations
  computed on both the initial merged TSS set (Stage=Initial) and the
  filtered final TSS set (Stage=Final), enabling comparison across filter
  stringencies.

### Changed
- `qc/qc_cs.txt` → `qc/qc_initial_cs.txt` and `qc/qc_in.txt` →
  `qc/qc_initial_in.txt` to make explicit that these metrics are computed
  on the initial (pre-filter) merged TSS set.
- `qc/stats_cs.txt` → `qc/stats_initial_cs.txt` and `qc/stats_in.txt` →
  `qc/stats_initial_in.txt` for the same reason.

### Fixed
- Sample sheet validation no longer fails when an optional column
  (`read_r2`, `input_name`, `replicate`) is present in the header but has
  all-empty values.  pandas reads empty cells as `NaN`; the `fillna`
  normalisation now runs before `snakemake.utils.validate` so the
  DataFrame is clean before the JSON schema type check.

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

[0.4.0]: https://github.com/bjmt/csRNAseq_pipeline/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/bjmt/csRNAseq_pipeline/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/bjmt/csRNAseq_pipeline/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/bjmt/csRNAseq_pipeline/releases/tag/v0.1.0
