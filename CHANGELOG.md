# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.16.2] - 2026-09-08

### Fixed
- `r-data-table` added to `workflow/envs/minimal_env.yaml` and `full_env.yaml`, and to
  `tsl/singularity.def`. `qc_complexity.R` has needed it since 0.16.0, so a fresh install
  from any of the shipped environment definitions could not run that rule. It passed on
  the development machine only because the host R library happened to carry the package.
- `bam2td` added to `tsl/singularity.def`, which had never installed it. `make_tagdir`
  calls it, so the container could not run the workflow past alignment. The build now
  fails at that point rather than at the first `make_tagdir` job if the binary does not
  land where expected.
- The compare links at the foot of this file pointed at the wrong GitHub organisation and
  stopped at 0.4.0. They now cover every release.

### Changed
- Documentation reorganised. The README is roughly a quarter of its former length and
  covers installation, running, the config keys worth knowing, the rules, the outputs and
  reproducibility. Three new files under `docs/` take the detail: `configuration.md` for
  the keys that need more than a schema description, `qc-metrics.md` for every QC column
  and what the numbers do not mean, and `why-these-defaults.md` for the reasoning behind
  the defaults and the gates.
- Neither the documentation nor this file refers to any particular dataset any more, and
  no example uses a real sample name.
- `tsl/singularity.def` also installs `gawk` explicitly, and `pyyaml` and `pytest` so the
  test suite can run inside the container.

### Note
- `workflow/Snakefile` declares `conda: "workflow/envs/env.yaml"`, but that file is listed
  in `.gitignore` and so is not in the repository. A `--use-conda` run from a fresh clone
  therefore has no environment to build. Either commit a portable version of that file or
  point the directive at `minimal_env.yaml`.

## [0.16.1] - 2026-09-07

### Changed
- Refreshed `tests/data/golden/qc_final_cs.txt` for the three columns 0.16.0 added,
  taking that table to 47 columns. No other value moved, and `qc_final_in.txt` is
  unchanged at 32 columns because the complexity metrics are csRNA-only. Two runs into
  separate directories agreed on every `*.complexity.txt` file and on both QC tables, so
  `qc_complexity` keeps the reproducibility property 0.15.0 established.

## [0.16.0] - 2026-09-07

### Added
- `TSSDistinctFrag`, `TSSChao1` and `TSSSaturation` in `qc_final_cs.txt`, from a new
  per-library `qc_complexity` rule. The QC table had no measure of library complexity,
  which is the binding constraint at low RNA input and which no other column catches.

  A fragment is a `(position, strand, length)` triple from the tag directory rather than
  a 5' position, because one in-TSS position carries several distinct lengths and a count
  of positions alone would treat much of that diversity as duplication. The count is
  restricted to the final TSS set, without which a library made of scattered background
  scores as more complex than a good one. Chao1 needs only singleton and doubleton
  counts, so it avoids the depth matching that would otherwise discard most of a good
  library's data, and Good's coverage answers whether more sequencing would pay. None of
  the three is gated on.
- `qc_five_prime` and `qc_complexity` in the validate layer's expected rules. The former
  had been missing since 0.11.0.

## [0.15.0] - 2026-09-04

### Changed
- Refreshed `tests/data/golden/qc_final_*.txt` for the read-length columns 0.13.0 and
  0.14.0 added, taking the csRNA table to 44 columns and the input table to 32. Every
  shared value is unchanged apart from `MinReplCorrPearson`, whose last-bit drift the
  ordering fix below settles.

### Fixed
- Two runs on the same data now write the same bytes. Every table HOMER writes came out
  in the order its threads finished, and nothing downstream imposed one, so a rerun
  reshuffled the rows of otherwise identical results. The per-library TSS tables and the
  BED derived from them, the four initial quantification tables, both final
  quantifications and the two tables `normalize_tss_quantification` writes are now sorted
  on a coordinate key by a shared `workflow/scripts/sort_table.sh`, under `LC_ALL=C` with
  a tab separator so the order follows neither the host's locale nor `sort`'s whitespace
  splitting.

  Nothing read those tables positionally, so no result depended on the order. Two things
  did follow from it: the count tables disagreed with `tss.final.bed`, which keeps the
  consensus BED's order throughout, and `cor()` added the same pairs up in a different
  order each run, which moved the Pearson replicate correlation in the last bit of a
  double. That is what had been changing `MinReplCorrPearson` on every golden refresh
  since 0.11.0.
- `makeUCSCfile` chose a random RGB triple per run, so all ten raw bedGraphs differed on
  their `track` line while their data lines matched. The rule now passes `-color`.

### Note
- What still differs between two runs is timing and provenance rather than results: the
  `*.benchmark.txt` records, the tool logs, `run_info.txt`, and one line of
  `tss/*.stats.txt` where HOMER records its own random draw. Run into a different
  directory and the files that record a tool's own command line differ too.

## [0.14.0] - 2026-09-04

### Added
- `ModeReadFraction`, the percentage of reads at `ModeReadLength`, so the height of the
  tallest bin rather than which bin it is. Of the read-length columns this is the one
  that separates library types, because a library made of one processed species
  concentrates in a few lengths while initiation over many promoters does not. Reported
  rather than gated on, since it largely restates `PctNuclear` and `csFRiP`.

### Changed
- Corrected the 0.13.0 claim for `P20ReadLength` and `P80ReadLength`. Neither the
  percentiles nor the spread between them separates csRNA libraries from input libraries.
  What holds is the paired comparison of a csRNA library against its own matched input.
  The columns stay, described as a summary of the distribution's shape rather than a test.

## [0.13.0] - 2026-09-04

### Added
- `P20ReadLength` and `P80ReadLength` beside `MedianReadLength`, from the same HOMER
  histogram and the same pass, so the shape of a library's length distribution is
  readable and not just its centre. The 80th percentile alone does not do the job,
  because size selection moves every percentile together. See 0.14.0, which corrects what
  this entry originally claimed for the pair.

## [0.12.2] - 2026-09-04

### Changed
- Refreshed `tests/data/golden/qc_final_*.txt`, which 0.12.0 flagged as no longer
  matching the column set. They gain `MedianReadLength`, `ModeReadLength`, `snRNA5pPct`
  and `tRNA5pPct` and lose `StatusBasis`. Every shared value is unchanged except
  `MinReplCorrPearson` in the last bit of a double, so this is a column change rather
  than drift in any metric.

## [0.12.1] - 2026-09-03

### Changed
- Corrected what 0.12.0 claimed about `ModeReadLength`. A healthy csRNA library is nearly
  flat in length, so the winning bin flips between replicates of one condition on
  differences too small to mean anything. The column still earns its place on libraries
  that are not flat, where a discrete population does take the mode, but what separates
  those cases is the height of the top bin, which this column does not carry.
  `MedianReadLength` is the stable summary of the two.

## [0.12.0] - 2026-09-03

### Added
- `MedianReadLength` and `ModeReadLength` in all four QC tables, input libraries
  included. Both come from HOMER's `tagLengthDistribution.txt`, so they describe the
  reads that survived filtering into the tag directory rather than everything the trimmer
  emitted, and they cost no extra pass over the data. The average HOMER prints in that
  file's header is pulled toward the tail by a few long reads.
- README rows for `snRNA5pPct` and `tRNA5pPct`, which 0.11.0 added without documenting.

### Removed
- `StatusBasis`. It repeated one constant on every row to say which table you were
  reading, which the filename already says. What it was carrying, that
  `collect_consensus_tss` reads the initial `Status` and nothing reads the final one, is
  now a comment beside each `Status` assignment and a paragraph in the documentation.

## [0.11.0] - 2026-09-03

### Added
- `qc/snrnas`, a path to a stranded BED of snRNA loci, and two columns in
  `qc_final_cs.txt`: `snRNA5pPct` and `tRNA5pPct`, the share of signal within 500 bp of
  an annotated 5' end that sits within 5 bp of it. A new per-library `qc_five_prime` rule
  computes them from the raw bedGraphs, whose coordinates are already read 5' ends.

  Sm-class snRNAs are capped Pol II transcripts, so a working library puts almost all of
  their signal on the 5' end, while tRNA 5' ends carry a monophosphate from RNase P and
  become ligation-competent as soon as the phosphatase step fails. `snRNA5pPct` is the
  one to gate on; `tRNA5pPct` moves only in the most extreme cases and is reported as a
  cross-check on which chemistry failed. Combining the two into a ratio or a difference
  is worse than `snRNA5pPct` alone.

  Each is a ratio computed inside one library against its own neighbourhood, so unlike
  `PhosEfficiency` they need no matched input and no normalisation. `qc/trnas` is reused
  for the tRNA half and now wants a strand column; a BED without one is skipped with a
  warning rather than scored wrongly. Off by default: leave `qc/snrnas` blank and only
  `tRNA5pPct` is computed.

## [0.10.1] - 2026-08-27

### Added
- `program/keep_below_mapq_sample`, an integer, default 0: keep this many reads per
  library that aligned but fell below `filtering/alignment_mapq`, as
  `qc/<sample>.belowmapq.sample.fastq.gz`. `keep_unmapped_sample` samples only records
  carrying the unmapped flag, so the multi-mapping fraction, which is usually the larger
  class, was counted and then discarded unseen. Separate key and file from the unmapped
  sample, because a read that never aligned and a read that placed equally well elsewhere
  mean opposite things. Sequences are written as sequenced, reverse-complemented back
  where the alignment was on the reverse strand.

### Fixed
- `merge_initial_tss` declared its per-sample BED inputs as a bare generator expression
  rather than a list. A generator yields its items once, so the rule's dependency set came
  back one library short, always the first in the sample table. The shell builds its file
  list from `params`, so the merge never omitted a library; what it lost was the guarantee
  that the library had finished being called. A rerun over existing BED files is the case
  to worry about, since `cat` can then read a file that is being rewritten. Present since
  the first Snakemake version.
- `keep_unmapped_sample` wrote a 0-byte and therefore unreadable `.fastq.gz` for a library
  with no unmapped reads, because the gzip pipe closed whether or not anything had been
  written to it. A class with no reads now leaves no file, and the zero is recorded in
  `qc/<sample>.aln.raw.txt` as before.

## [0.10.0] - 2026-08-26

### Added
- `qc/min_log2_fold`: a library whose TSS calling used an enrichment threshold below this
  fails QC. Unset by default, falling back to `program/homer/tss/default_log2_fold`.
  With an annotation HOMER chooses the threshold per library from the data, writes it to
  `tss/<sample>.stats.txt`, and offers no flag to floor it. For a library that has lost
  its capped signal it can come out at or below zero, which accepts TSSs carrying no more
  signal than their own input, and because the consensus is a union those calls propagate
  downstream. The pipeline reads the threshold back and gates on it, reporting it as a new
  `Log2FoldThreshold` column in both QC tables.
- `AlignedReads`, `BelowMapqReads` and `FilteredOutReads` in the stats and QC tables,
  between `TrimmedReads` and `TotalReads`. That gap is the largest single loss anywhere in
  the pipeline and nothing recorded where it went. Three things could account for it with
  opposite implications: reads that never aligned, reads below `filtering/alignment_mapq`,
  and reads dropped by `exclude_flags`, `min_alignment_length` or `max_mismatch`. The
  counts come from a pass-through awk in the alignment pipe, since the unfiltered
  alignment never reaches disk and the flagstat the pipeline keeps runs on the filtered
  BAM, where it reports 100% mapped for every library. Per-library detail, including a
  MAPQ histogram, lands in `qc/<sample>.aln.raw.txt`.
- `program/keep_unmapped_sample`, an integer, default 0: keep that many unmapped reads per
  library as `qc/<sample>.unmapped.sample.fastq.gz`. With `keep_bam` and
  `keep_trimmed_fastq` off, nothing about the unaligned fraction survives a run, so a
  contaminant screen afterwards means aligning again.
- STAR's `Log.final.out` is kept rather than deleted, appended to
  `qc/<sample>.aligner.log`. STAR keeps unmapped reads out of its BAM entirely, so its own
  report is the only place its read accounting exists.
- A documented caveat that every per-class QC figure (`PretRNAPct`, `miRNAPct`, any rRNA
  number) comes from the filtered BAM and therefore measures the unique-mapping tail of
  those species rather than their abundance.
- `qc/min_cs_frip_final` and `qc/min_pct_nuclear_final`, gates for the final QC table.
  Both tables report `csFRiP` but measure it against different regions, so one number
  cannot serve both: every library's final `csFRiP` falls when the consensus shrinks, so a
  gate tuned on the initial table goes on to flag sound libraries in the final one. Unset
  by default, falling back to `min_cs_frip` and `min_pct_nuclear`.
- `StatusReason` in both csRNA QC tables, naming the gates that tripped. `Status` alone
  said only Ok or FAIL, which is opaque once more than one gate is in play.
- `StatusBasis` in both csRNA QC tables, recording that `collect_consensus_tss` consumes
  the initial `Status` and nothing consumes the final one. Removed again in 0.12.0.

### Changed
- **`qc/min_log2_fold` changes which libraries pass QC**, and with
  `filtering/exclude_failed_from_consensus` on it changes the consensus set, so an
  existing analysis will shift. Compare before adopting, as with the `-gtf` entry in
  0.8.0.
- `rule align` runs a counter in the alignment pipe, so its shell command has changed and
  Snakemake will want to realign existing runs. The alignment itself is untouched: the awk
  passes the stream through byte for byte, which the test suite asserts.
- `find_tss_initial` post-processes `tss/<sample>.stats.txt`, which is now a declared
  output of that rule rather than an untracked side effect. HOMER writes the
  stable-transcript and promoter-distal fractions whether or not the input they need was
  supplied, so both read as measurements; each is rewritten to `na` when HOMER's own report
  shows it had nothing to compute it from. A real measurement is never overwritten, and
  anything parsing that file should expect `na` in those two positions.
- `collect_consensus_tss.R` sends its warnings to stderr rather than stdout, so the notice
  that a failed library still contributes to the union is no longer lost among the
  progress narration in the Snakemake log.

### Fixed
- Editing `qc/min_cs_frip` or `qc/min_pct_nuclear` and rerunning was a silent no-op. Both
  reached the QC scripts through `snakemake@config`, which no rerun trigger can see, so
  the stale `Status` stayed on disk with nothing in the log to say the new gate had been
  ignored. Both now arrive as rule params, which `--rerun-triggers params` acts on.
  `collect_consensus_tss.R` had the same bug for `filtering/tss_min_reps`. No
  `snakemake@config` read is left in any script.
- A QC gate that cannot be evaluated, an `NA` `csFRiP` from a library with no nuclear
  reads, left `Status` as `NA`, which `collect_consensus_tss` treated as neither Ok nor
  FAIL. It now counts as a failure, named in `StatusReason`.
- The documentation pointed at `qc/qc_final_cs.txt` and `qc/qc_final_in.txt`; both have
  been at the output root since the output paths were reorganised.

## [0.9.0] - 2026-08-04

### Removed
- `program/homer/tss/rnaseq_tagdir` and `rnaseq_stability_only`, so HOMER's `-rna` is no
  longer wired in. `-gtf` stays and carries the whole benefit anyway: on the test genome,
  `-gtf` alone yields the same true and false positive sets (15 and 5) and the same chosen
  input threshold (0.415) as `-gtf` with `-rna`, because the RNA-seq side only adds the
  stable-transcript filter on top. Dropping it removes an option that required a
  hand-built HOMER tag directory from a splice-aware aligner, which this pipeline's trim
  and alignment settings cannot produce, along with a strand-orientation trap that would
  have made the stability call meaningless if got wrong.
- The end-to-end step and wiring assertions now cover `-gtf` only.

### Fixed
- The documentation no longer cites `Skipping TSS assignment (can't find file for genome
  ...)` as evidence that the annotation is missing. That line comes from the
  `annotatePeaks` calls inside HOMER wanting a full genome directory and appears whether
  or not `-gtf` was given. The usable signal is the TP and FP counts.

## [0.8.1] - 2026-08-04

### Added
- `program/homer/tss/rnaseq_stability_only` passes HOMER's `-noFilterRNA`, so RNA-seq can
  fill the stable-transcript columns without filtering TSSs. Supplying `rnaseq_tagdir`
  does filter by default, which is HOMER's behaviour and was not documented. Default
  False, matching HOMER.
- End-to-end coverage for both optional annotation inputs, plus a GTF fixture for the
  synthetic genome. `tests/e2e/check_annotation_wiring.py` asserts on HOMER's own report
  rather than on the flags being present, since both options are easy to wire up so that
  they look fine and do nothing. It checks that HOMER read the GTF, built non-empty TP and
  FP sets from it (15 and 5 on the fixture), chose the input threshold from the data
  rather than falling back to `-defaultLog2Fold`, and applied a real RNA-seq threshold.

### Changed
- Documented that `rnaseq_tagdir` is a tag directory rather than FASTQ, and why the
  pipeline cannot build it: RNA-seq reads are spliced and longer, while the trim step caps
  reads at `max_read_length` and the aligners are configured for short unspliced reads.

### Fixed
- `program/homer/tss/rnaseq_tagdir` and `gtf` were accepted without any check beyond
  existence, so a path that was not a HOMER tag directory, or not a GTF, passed start-up
  and then failed deep inside HOMER. Both are validated up front now, the GTF for nine
  tab-separated fields and the tag directory for a `tagInfo.txt`, with errors that say
  what the option wants.

## [0.8.0] - 2026-08-04

### Added
- `program/homer/tss/gtf`: optional gene annotation passed to HOMER as `-gtf`. Without it
  HOMER cannot build the annotated-TSS and exon sets it uses as true and false positives
  to choose each library's enrichment threshold, so it falls back to `-defaultLog2Fold`
  and writes the promoter-distal and stable-transcript columns as placeholders. Empty by
  default, because the pipeline must keep working without an annotation, but it is the
  largest single quality lever in the configuration. Setting it changes TSS calls, so
  compare against an existing consensus before adopting it.
- `program/homer/tss/program`: choice of HOMER caller. `findcsRNATSS.pl` announces itself
  as a legacy placeholder for `findcsRNATSR.pl`; the two are the same program, their
  diff being renames and help text, and a row-for-row check on a full-size library found
  their tables identical. The pipeline normalises TSR's `.tsr.txt` and `.alltsr.txt` names
  to the `.tss.txt` the rest of it expects. The default stays on the legacy name because
  TSR ships an unpatched `use lib` pointing at its author's install, so it resolves
  `HomerConfig` only when HOMER's bin is on `PERL5LIB`.
- `program/homer/tss/pseudo_count`, `default_log2_fold` and `local_fold` expose HOMER's
  `-pseudoCount`, `-defaultLog2Fold` and `-L`, all defaulting to HOMER's own values.
  `pseudo_count` is what lets a TSS with no input coverage pass the enrichment test
  automatically. `-cpu` is now passed so HOMER uses the threads the rule reserves.
- `filtering/exclude_failed_from_consensus`, default False: drop csRNA libraries whose
  initial QC `Status` is FAIL from the consensus union. The union is a union, so a
  contaminated library's spurious TSSs otherwise enter the shared set and everything
  downstream quantifies them. A warning naming the libraries is printed either way.
- `run_info.txt`: per-run provenance. Pipeline version and git commit, flagged when the
  working tree is dirty, plus host, resolved config, sample-table checksum, and the
  version of every tool the run used.
- Golden-value regression test. `tests/data/golden/qc_final_*.txt` are compared
  value-by-value after the default e2e run, so a metric changing silently fails the suite,
  which the structural assertions never noticed. Refresh deliberately with
  `tests/scripts/update_golden.sh <e2e outdir>`. Skipped for the per-aligner runs, where a
  different aligner could legitimately shift values.

### Changed
- `tss_size_composition` is split into a per-library rule and a join, so Snakemake runs
  the tag-directory pass for each library in parallel instead of looping over them in one
  job.
- `tss.consensus.sizes.txt` carries `.top1` to `.topN` up to the new
  `filtering/tss_top_sizes_max`, default 5, instead of a single `.topn`, so
  `tss_top_sizes_n` can be retuned anywhere in that range without re-reading the tag
  directories. Keeping the full per-TSS length histogram would serve any n exactly but
  runs to millions of rows per library, against a handful of columns here.
- Documented what the QC table does not mean: that `csFRiP` and its neighbours measure
  different regions in the initial and final tables, that the contamination percentages
  are signal-relative rather than library-relative, that a high duplicate rate can mean a
  *better* library because concentrated reads duplicate more, and that `csEnrichment` does
  not discriminate at two to three libraries per condition.

## [0.7.0] - 2026-08-04

### Added
- Top-lengths filter, a second read-size composition filter that needs no prior knowledge
  of which lengths are contaminating. `tss_max_top_sizes_fraction`, default 1 and so
  disabled, drops a TSS when more than that fraction of its reads fall in its
  `tss_top_sizes_n` commonest read lengths, 2 by default. Genuine initiation is
  heterogeneous and spreads over tens of lengths, while a discretely processed RNA with
  fixed 5' and 3' ends puts nearly everything into one or two. It catches contaminants
  that fall outside any small RNA size class. Gated by the existing `tss_srna_min_reads`
  and `tss_srna_min_samples`, whose descriptions now say they govern both size filters.
- `tss.consensus.sizes.txt` gains a `.topn` column per library. Columns are written on
  demand: `.reads` always, `.srna` only when `tss_srna_sizes` is set, `.topn` only when
  the top-lengths filter is on.
- Unit coverage for the new filter: two more `normalize` invocations and a third
  `size_composition` invocation checking the `.topn` column layout. The tag fixtures now
  give most clusters a spread of read lengths, so both filters have real negatives as well
  as positives.

### Changed
- `qc_final_cs.txt` and `qc_final_in.txt` no longer carry `FinalRiP`, `FinalFRiP` and
  `FinalEnrichment`, which were literal copies of `csRiP`, `csFRiP` and `csEnrichment` in
  the same table. `FinalTSSDetected` stays, being a distinct metric. The duplication was
  hiding a naming subtlety, now documented: in the final tables those three names are
  computed against `tss.final.bed`, and in the initial tables against the initial merged
  set, so the two files must not be compared column by column.
- `miRNAPct` uses the same denominator as `PretRNAPct`, contaminant reads over contaminant
  plus reads-in-TSS, instead of contaminant over nuclear reads. The two contamination
  percentages sit side by side and were on different bases. Values shift up slightly;
  ratios between libraries are essentially unchanged.
- `tss_srna_min_reads` default raised from 30 to 100. A handful of reads gives a poor
  estimate of size composition, so the low floor cost specificity at weakly expressed
  TSSs. Genuine single-locus contaminants are abundant, so the higher floor keeps them.
  Only affects runs that enable a size filter; `tests/data/config.yaml` pins 30 so the
  unit fixtures still exercise the floor.
- Corrected the documented rationale for the small-RNA size filter. It previously said
  enrichment over the input fails to reject siRNA TSSs where the input is too shallow to
  measure a background. Shallow input turns out to be a minor sub-mode. What the filter
  really guards against is loss of library quality, since a library that has lost its
  capped signal is proportionally richer in siRNA than its own input, while read size does
  not degrade that way. Documentation only, no behaviour change.

### Fixed
- `gather_stats.sh` summed organelle tag counts in bash integer arithmetic. HOMER writes
  tag counts as floats, and any fractional total comes back from awk in `%.6g` scientific
  notation, which bash rejects outright, failing the rule. The sum and rounding now happen
  inside awk.
- `collect_consensus_tss.R` could place a TSS outside its chromosome, because widening to
  150 bp and the overlap-resolution shifts were unclamped and no seqlengths were set. Such
  a TSS is now slid back inside, preserving its width, and one wider than its own
  chromosome is left alone and reported.
- `collect_consensus_tss.R` picked one of each symmetric overlap pair with
  `seq(1, length(hits), by = 2)`, which assumes the two members land adjacent in the
  `Hits` object. That holds for equal-width sorted ranges but is not guaranteed, so it now
  selects on `queryHits < subjectHits`.
- `tss_size_composition.R` built its per-TSS, per-length key in integer arithmetic, which
  overflows silently to `NA` for a large genome sequenced with long reads. Now double.
- A failed `findcsRNATSS.pl` call is no longer silent. The rule still writes an empty TSS
  set so one bad library cannot stall a run, but it says so on stderr and in the log,
  because the only other symptom was the library quietly vanishing from the consensus
  union and a confusing `bedtools` error downstream.

## [0.5.0] - 2026-08-04

### Added
- Read-size composition filter on the consensus TSS set, driven by four new `filtering`
  keys: `tss_srna_sizes` (an array of read lengths), `tss_max_srna_fraction` (0.5),
  `tss_srna_min_reads` (30) and `tss_srna_min_samples` (1). A TSS is dropped from
  `tss.final.bed` when more than `tss_max_srna_fraction` of its reads fall in the listed
  sizes, in at least `tss_srna_min_samples` csRNA libraries that have at least
  `tss_srna_min_reads` reads there. Defaults to `[]`, which disables the filter.

  This catches uncapped small RNAs that the other filters cannot reach. In plants the
  21-25 nt siRNAs, above all the 24 nt Pol IV class over transposons, survive TEX and AP
  well enough to be called as TSSs. Read length separates them from genuine initiation,
  which is not confined to one small size class.
- New rule `tss_size_composition` and script `workflow/scripts/tss_size_composition.R`,
  measuring per-TSS read counts and small-RNA-sized read counts from the tag directories
  into `tss.consensus.sizes.txt`. The length column of the tag files supplies the sizes,
  so no extra HOMER run is needed. When `tss_srna_sizes` is empty the rule writes the TSS
  list and does not read the tag directories at all.
- Unit test `tests/unit/size_composition.smk` with tag-directory fixtures, plus two
  further `normalize` invocations covering the new filter at `tss_srna_min_samples` 1 and
  2.

### Changed
- `normalize_tss_quantification` stops with an explanatory message when the filters leave
  no TSSs at all, rather than failing inside edgeR with "'counts' must contain at least
  one value".

### Fixed
- `tests/scripts/generate_unit_fixtures.py` writes the `RawReads` and `TrimmedReads`
  columns added in 0.3.0, so the `qc_initial` and `qc_final` unit tests pass again.
- `tests/run_tests.sh` propagates unit-test failures out of `run_unit_smk` instead of
  returning the status of its last conditional, which had been reporting failed assertions
  as passes.

## [0.4.0] - 2026-05-13

### Changed
- `make_tagdir` uses `bam2td` instead of HOMER `makeTagDirectory`. bam2td reads sorted
  BAMs directly, with no intermediate SAM conversion, and produces a compatible HOMER tag
  directory, so `shadow: "minimal"` and the `samtools view -h` step are gone.
  `makeTagDirectory -update -checkGC -genome` runs afterwards to populate the QC files
  bam2td does not generate. The `homer.tagdir.extra_args` key is still supported and goes
  to bam2td.

## [0.3.0] - 2026-05-13

### Added
- `RawReads` and `TrimmedReads` columns in `stats_initial_cs.txt` and
  `stats_initial_in.txt`, propagated to all QC output files, parsed from bfqutils
  `.trimming.txt` logs. Multi-file samples are handled by summing across their log
  entries.
- `miRNAPct` column, miRNA reads over `NuclearReads` times 100, in all QC output files,
  computed alongside the existing `miRNA` and `miRNADepletion` metrics.
- `tss_min_cs_in_ratio` and `tss_min_ratio_samples`: a csRNA/input read-ratio filter
  applied in `normalize_tss_quantification` after the CPM filter. It uses the
  already-computed input quantification of the initial csRNA TSS set, so no extra HOMER
  run is required. Defaults to 0, disabled.

### Changed
- The `key` column is renamed `group` in all QC output files.
- `qc_final_cs.txt` and `qc_final_in.txt` are written to the results root rather than
  `qc/`.
- The `tss_merged/` outputs and all files derived from them land in the results directory
  rather than the intermediate directory.
- `generate_normalized_bw` memory raised from 2,000 MB to 4,000 MB.

## [0.2.0] - 2026-05-06

### Added
- End-to-end tests for bwa-mem and hisat2, runnable as `./tests/run_tests.sh bwa-mem` and
  `./tests/run_tests.sh hisat2`, and included in the `all` subcommand.
- `workflow/envs/minimal_env.yaml` rewritten as a portable, cross-platform conda spec with
  the correct channel order, minimum-version constraints, and notes for the manual-install
  dependencies.
- `workflow/envs/full_env.yaml` for development and the full test suite, adding STAR,
  bowtie2, hisat2, pytest, snakefmt and vim on top of the minimal spec.
- `qc_final_cs.txt` and `qc_final_in.txt` carry the full set of QC metrics, re-derived
  from the raw quantification files rather than copied from the initial QC step, so they
  reflect the filtered TSS set.
- `qc/replicate_correlation.txt` carries correlations computed on both the initial merged
  TSS set (`Stage=Initial`) and the filtered set (`Stage=Final`), so filter stringencies
  can be compared.

### Changed
- `qc/qc_cs.txt` and `qc/qc_in.txt` become `qc/qc_initial_cs.txt` and
  `qc/qc_initial_in.txt`, and `qc/stats_*.txt` likewise, to make explicit that those
  metrics are computed on the initial pre-filter merged TSS set.

### Fixed
- Sample sheet validation no longer fails when an optional column (`read_r2`,
  `input_name`, `replicate`) is present in the header with all-empty values. pandas reads
  empty cells as `NaN`, so the `fillna` normalisation now runs before
  `snakemake.utils.validate`.

## [0.1.0] - 2026-05-04

### Added
- Snakemake workflow with 13 rules covering the full analysis: adapter trimming
  (bfqutils), genome alignment (bwa-aln, bwa-mem, STAR, bowtie2, hisat2), HOMER tag
  directory creation, initial TSS detection, strand-separated bedGraph generation,
  alignment statistics, TSS merging, quantification, TMM normalisation, consensus TSS
  collection, CPM filtering, post-filter QC, and normalised bigWig generation.
- QC metrics: FRiP, nuclear read fraction, csRNA enrichment over input, miRNA and tRNA
  contamination, Pearson replicate correlation, and per-sample alignment statistics.
- `input_name` column in the sample table, so several conditions can share one input
  library.
- CPM filter (`tss_min_cpm`, `tss_min_samples`) to remove low-confidence TSSs from the
  consensus set, with a dedicated post-filter QC step.
- Small-RNA masking option (`mask_srnas_from_bw`) to exclude miRNA and tRNA loci from the
  normalised bigWig tracks.
- Schema validation for the sample table and configuration file via
  `snakemake.utils.validate`.
- `conda` environment definition (`workflow/envs/env.yaml`).
- Test suite under `tests/` with four layers: schema validation, Snakemake dry-run
  validation, R-script unit tests, and end-to-end tests against a synthetic 80 kb genome.

### Fixed
- STAR alignment: bfqutils appends a trailing space to every read ID line, which STAR
  2.7.x rejects. Reads are piped through `sed` inside a bash process substitution to strip
  trailing whitespace before STAR sees them.
- STAR output prefix derivation: `Path.with_suffix()` on a Snakemake placeholder string
  corrupted the closing brace. The prefix is computed at shell runtime with `dirname` and
  `basename`.
- `generate_normalized_bw.R`: negative-strand bedGraphs are no longer empty when every
  read is forward-strand, which had crashed `rtracklayer::export.bw()`.
- `zcat` on macOS cannot decompress `.gz` files; replaced with `gzip -dc` throughout.

[Unreleased]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.16.2...HEAD
[0.16.2]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.16.1...v0.16.2
[0.16.1]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.16.0...v0.16.1
[0.16.0]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.15.0...v0.16.0
[0.15.0]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.14.0...v0.15.0
[0.14.0]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.13.0...v0.14.0
[0.13.0]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.12.2...v0.13.0
[0.12.2]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.12.1...v0.12.2
[0.12.1]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.12.0...v0.12.1
[0.12.0]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.11.0...v0.12.0
[0.11.0]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.10.1...v0.11.0
[0.10.1]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.10.0...v0.10.1
[0.10.0]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.9.0...v0.10.0
[0.9.0]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.8.1...v0.9.0
[0.8.1]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.8.0...v0.8.1
[0.8.0]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.5.0...v0.7.0
[0.5.0]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/noborilab/csRNAseq_pipeline/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/noborilab/csRNAseq_pipeline/releases/tag/v0.1.0
