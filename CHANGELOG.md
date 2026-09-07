# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.16.1] - 2026-09-07

### Changed
- Refreshed `tests/data/golden/qc_final_cs.txt` for the three columns 0.16.0 added,
  `TSSDistinctFrag`, `TSSChao1` and `TSSSaturation`, which takes that table to 47
  columns. Nothing else in it moved. `qc_final_in.txt` is byte-identical and stays at 32,
  since the complexity metrics are csRNA-only.

  Two runs into separate output directories agreed on all three `*.complexity.txt` files
  and on both QC tables, so the new `qc_complexity` rule holds the reproducibility
  property 0.15.0 established, and `MinReplCorrPearson` sat on 0.5273890876768268 in both
  rather than drifting a bit per run as it did before that release.

## [0.16.0] - 2026-09-07

### Added
- `TSSDistinctFrag`, `TSSChao1` and `TSSSaturation` in `qc_final_cs.txt`, from a new
  per-library `qc_complexity` rule. The QC table has had no measure of library
  complexity, which on the protocol-development panel is the binding constraint at low
  RNA input: the 1 ug libraries have the highest csFRiP in that panel and a third of the
  distinct molecules, so nothing already in the table would have caught them.

  The unit is the `(position, strand, length)` triple the tag directories already carry,
  not the 5' position the bedGraphs carry. That is not a detail: measured on this panel,
  5' position alone treats between a half and three quarters of the distinct fragments as
  duplicates, because an in-TSS position holds a median of about 3.7 distinct lengths.

  Restricted to the final TSS set, which is load-bearing rather than cosmetic. Measured
  genome-wide, the two libraries that failed QC outright ranked first and second on
  complexity, since a library made of scattered background has more distinct positions
  than a good one. Restricted to called TSSs they fall to ninth and eighteenth of
  twenty-two.

  Chao1 rather than a fixed-depth subsample because the shallowest library on the panel
  had 2.5M in-TSS reads: matching depth would either throw away 90% of a good library's
  data or return `NA` for exactly the libraries worth judging. It agrees with a
  depth-matched count at rho = 0.88 among libraries that pass the other gates.

  None of the three is gated on. `TSSChao1` still inflates for background-heavy
  libraries, and `TSSSaturation` runs at rho = 0.77 against `csFRiP` among passing
  libraries, so it partly restates purity. They are reported to be read beside it.

### Note
- `tests/data/golden/qc_final_*.txt` gain three columns and need
  `tests/scripts/update_golden.sh <e2e outdir>` on the next e2e run.

## [0.15.0] - 2026-09-04

### Changed
- Refreshed `tests/data/golden/qc_final_*.txt` for the read-length columns 0.13.0 and
  0.14.0 added, `P20ReadLength`, `P80ReadLength` and `ModeReadFraction`, which takes the
  csRNA table to 44 columns and the input table to 32. Every value in a column shared
  with the previous copy is unchanged, apart from `MinReplCorrPearson`, whose last-bit
  drift the row-ordering fix below settles.

### Fixed
- Two runs of the pipeline on the same data now write the same bytes. Every table HOMER
  produces was coming out in the order its threads happened to finish in, and nothing
  downstream imposed one, so a rerun reshuffled the rows of results that were otherwise
  identical. Sorted on the way out of the rule that makes each one, all on the same
  coordinate key, by a new `workflow/scripts/sort_table.sh`:

  | Table | Rule |
  |-------|------|
  | `tss/{sample}.tss.txt`, `tss/{sample}.alltss.txt` | `find_tss_initial` |
  | `bed/{sample}.tss.bed` | `find_tss_initial`, inherited through `pos2bed.pl` |
  | `tss_merged/all_cs.tss_merged_quant_{cs,in}.txt` | `quantify_initial_cs_tss` |
  | `tss_merged/all_in.tss_merged_quant_{cs,in}.txt` | `quantify_initial_in_tss` |
  | `tss.consensus.homer.raw.txt` | `quantify_final_tss` |
  | `tss.final.in.raw.txt` | `quantify_final_in_tss` |
  | `tss.final.raw.txt`, `tss.final.cpm.txt` | `normalize_tss_quantification`, ordered by `tss.final.bed` |

  Nothing read any of these positionally, so no result ever depended on the order, but two
  things followed from it. The count tables disagreed with `tss.final.bed`, which keeps the
  consensus BED's order throughout, and `cor()` added the same pairs up in a different
  order on every run, which moved both the `Initial` and the `Final` Pearson correlation
  in the last bit of a double. That is what has been quietly changing
  `MinReplCorrPearson` on every golden refresh for three releases. It is now stable at
  0.5273890876768268, and `MinReplCorrSpearman` never moved because rank sums are exact.

  `sort` runs under `LC_ALL=C` with a tab separator, so the order does not follow the
  machine's locale and lines up with the columns being named rather than with whitespace.

- `makeUCSCfile` was picking a random RGB triple for each raw bedGraph's `track` line, so
  all ten bedGraphs differed between runs on that one line while their 247 data lines were
  identical. The rule now passes `-color`, blue for the forward strand and red for the
  reverse.

### Note
- What is left is timing and provenance rather than results. Two runs into the same output
  directory now agree on 190 of 250 files, and the rest are the 46 `*.benchmark.txt`
  wall-clock records, 11 tool logs, `run_info.txt`, which is meant to be unique per run,
  and the three csRNA `tss/*.stats.txt`, which differ only on HOMER's own
  `Random number used:` line. That line is left as HOMER writes it. Run into a different
  directory and a further 29 files differ only where a tool records its own command line,
  absolute paths included.

## [0.14.0] - 2026-09-04

### Added
- `ModeReadFraction`, the percentage of reads sitting at `ModeReadLength`. Height of the
  tallest bin rather than which bin it is, from the same histogram and the same pass.

  This is the read-length column that carries information. Across the 34 protocol-
  development libraries every passing csRNA library sat between 4.8 and 8.6, every
  failing one between 9.9 and 15.5, and every input library between 10.9 and 26.1, with
  no overlap between any of the three. A library made of one processed species
  concentrates in a few lengths, and initiation spread over thousands of real promoters
  does not, which is the same reasoning behind `filtering / tss_max_top_sizes_fraction`
  applied to the library rather than to a cluster.

  It is reported, not gated. On that panel it runs at r = -0.92 against `PctNuclear` and
  -0.74 against `csFRiP`, so it mostly restates gates the table already has, and the
  margin between the passing libraries and the rest is 1.3 percentage points on 22
  csRNA libraries from one experiment.

### Changed
- Corrected the 0.13.0 claim for `P20ReadLength` and `P80ReadLength`. Measured on all 34
  tag directories rather than the 8 in the cached table, neither the percentiles nor the
  spread between them separates csRNA libraries from input libraries: P80 minus the
  median runs 5 to 12 in the csRNA libraries and 5 to 10 in the inputs. What does hold is
  the paired version, where a csRNA library carries the longer upper tail than its own
  matched input in 11 of 12 pairs, by 1 to 4 nt. The columns stay, described as what they
  are, which is a description of the distribution's shape rather than a test.

## [0.13.0] - 2026-09-04

### Added
- `P20ReadLength` and `P80ReadLength` beside `MedianReadLength`, from the same HOMER
  histogram and the same single pass, so the shape of a library's length distribution
  is readable and not just its centre.

  The 80th percentile alone does not do the job, because size selection moves every
  percentile together: on the protocol-development panel a 40-70 nt input library and a
  20-70 nt csRNA library both land on P80 = 39 while looking nothing alike.

  This entry originally claimed that P80 minus the median separated csRNA libraries from
  input libraries without overlap, on 8 libraries from a cached table. Checked against
  all 34 tag directories it does not. See 0.14.0.

## [0.12.2] - 2026-09-04

### Changed
- Refreshed `tests/data/golden/qc_final_*.txt`, which 0.12.0 flagged as no longer
  matching the column set. They now carry `MedianReadLength`, `ModeReadLength`,
  `snRNA5pPct` and `tRNA5pPct`, and no longer carry `StatusBasis`. Every value in a
  column the two copies share is unchanged, except `MinReplCorrPearson` in the last
  bit of a double, which is far inside the comparison's tolerance, so this is a column
  change and not a drift in any metric.

### Note
- `snRNA5pPct` and `tRNA5pPct` are NA for all three libraries in the refreshed golden
  tables, because `qc/snrnas` is unset in the test config and no synthetic read lands
  within 500 bp of the placeholder tRNA locus in `tests/data/trnas.bed`. The golden
  comparison therefore pins the two columns' presence but not their arithmetic. Giving
  them real values means putting signal at those loci, which would move `PretRNA`,
  `PretRNAPct` and `PhosEfficiency` as well.

## [0.12.1] - 2026-09-03

### Changed
- Corrected what 0.12.0 claimed about `ModeReadLength`. That entry said the mode landing
  on 24 nt showed a 24 nt siRNA population sitting under the initiation signal. Checking
  the underlying histograms does not support reading it that way. A healthy csRNA library
  is nearly flat in length, the top bin holding 5 to 6 per cent, so the winning length
  flips between replicates of the same arm (24 in one, 31 in the other) on differences
  too small to mean anything. The column is still worth reporting, because a library with
  a genuinely discrete population is not flat: an input library puts 22 per cent at 24 nt
  and a phosphatase-free arm 15 per cent at 30 nt, and there the mode does land on the
  population. What separates the two cases is the height of the top bin, which this column
  does not carry. `MedianReadLength` is the stable summary of the two.

## [0.12.0] - 2026-09-03

### Added
- `MedianReadLength` and `ModeReadLength` in all four QC tables, for input libraries as
  well as csRNA ones. Both come from HOMER's `tagLengthDistribution.txt`, so they
  describe the reads that survived filtering into the tag directory rather than
  everything the trimmer emitted, and they cost no extra pass over the data. The
  average HOMER prints in that file's header is pulled toward the tail by a few long
  reads. On the enzyme panel the median runs 24 to 36, input libraries sitting shorter
  than the csRNA libraries drawn from the same RNA. See 0.12.1 for what the mode does
  and does not support.
- README rows for `snRNA5pPct` and `tRNA5pPct`, which 0.11.0 added to the tables
  without documenting.

### Removed
- `StatusBasis`. It repeated one constant on every row of a table (`initial` or
  `final-advisory`) to say which table you were reading, which the table's own filename
  already says. The fact it was carrying, that `collect_consensus_tss` reads the initial
  `Status` and nothing reads the final one, is now a comment beside each `Status`
  assignment and a paragraph in the README, where it is not repeated 22 times per run.

### Note
- `tests/data/golden/qc_final_*.txt` do not match the current column set, and have not
  since 0.11.0 added `snRNA5pPct` and `tRNA5pPct`. Refresh them with
  `tests/scripts/update_golden.sh <e2e outdir>` on the next e2e run.

## [0.11.0] - 2026-09-03

### Added
- `qc/snrnas`, a path to a stranded BED of snRNA loci, and two columns in
  `qc_final_cs.txt`: `snRNA5pPct` and `tRNA5pPct`, the share of signal within 500 bp
  of an annotated 5' end that sits within 5 bp of it. Computed per library by the new
  `qc_five_prime` rule, which reads the raw bedGraphs, whose coordinates are already
  read 5' ends.

  Sm-class snRNAs are capped Pol II transcripts, so a working library puts almost all
  of their signal on the 5' end, while tRNA 5' ends carry a monophosphate from RNase P
  and become ligation-competent as soon as the phosphatase step fails.

  `snRNA5pPct` is the one to gate on. On the enzyme panel it separates the two groups
  without overlap: worst library with the phosphatase 92.4, best without it 88.9, a gap
  of 4.5 SD. `tRNA5pPct` does not separate on its own, spanning 30.3 to 65.9 among
  libraries that have the phosphatase, because it only moves in the most extreme arm;
  it is reported because it is free and says which chemistry failed, not because it
  gates. Combining the two into a ratio or a difference was tested and is worse than
  `snRNA5pPct` alone: both lose the clean split, because the tRNA noise swamps the
  snRNA signal.

  Both are ratios computed inside one library against its own neighbourhood, so
  unlike `PhosEfficiency` they need no matched input and no normalisation, and work on
  a sample that has none. `qc/trnas` is reused for the tRNA half; it now wants a strand
  column, and a BED without one is skipped with a warning rather than scored wrongly.
  Off by default: leave `qc/snrnas` blank and only `tRNA5pPct` is computed.

## [0.10.1] - 2026-08-27

### Added
- `program/keep_below_mapq_sample`, an integer, default 0: keep this many reads per library
  that aligned but fell below `filtering/alignment_mapq`, as
  `qc/<sample>.belowmapq.sample.fastq.gz`. `keep_unmapped_sample` samples only records
  carrying the unmapped flag, so the multi-mapping fraction, which is the larger class and
  the one that bears on where a big alignment loss went, was counted and then discarded
  unseen. Separate key and separate file from the unmapped sample, because a read that
  never aligned and a read that placed equally well elsewhere mean opposite things. Works
  for every aligner including STAR, whose BAM does contain its below-MAPQ reads.
  Sequences are written as sequenced, reverse-complemented back where the alignment was on
  the reverse strand.

### Fixed
- `merge_initial_tss` declared its per-sample BED inputs as a bare generator expression
  rather than a list. A generator yields its items once, so the rule's dependency set came
  back one library short, always the first in the sample table. The shell builds its file
  list from `params` rather than from `input`, so the merge never silently omitted a
  library; what it lost was the guarantee that the library had finished being called.
  Seen both ways: a hard failure where the BED did not exist yet, and, on a 34-library
  panel, the merge starting in the same second the last TSS call finished. A rerun over
  existing BED files is the case to worry about, since `cat` can then read a file that is
  being rewritten. Present since the first Snakemake version. Every other rule builds its
  file lists with `list()`; this one now does too.
- `keep_unmapped_sample` wrote a 0-byte, and therefore unreadable, `.fastq.gz` for a
  library with no unmapped reads, because the gzip pipe was closed whether or not anything
  had been written to it. A class with no reads now leaves no file, and the zero is
  recorded in `qc/<sample>.aln.raw.txt` as before.

## [0.10.0] - 2026-08-26

### Added
- `qc/min_log2_fold`: a library whose TSS calling used an enrichment threshold below this
  now fails QC. Unset by default, falling back to `program/homer/tss/default_log2_fold`,
  so with the shipped default of 1 any library whose threshold came out below 2-fold is
  flagged; see the warning under Changed, because this one moves TSS calls.
  With an annotation HOMER chooses the threshold per library from the data, writes it to
  `tss/<sample>.stats.txt`, and offers no flag to floor it. On libraries that have lost
  their capped signal it picks values at or below zero (0.144, 0.007, -0.056 and -0.680 on
  one 22-library panel), which accepts clusters carrying less signal than their own input.
  Their promoter-distal fractions ran 48 to 54% against 19 to 24% for the sound libraries,
  so the extra calls are mostly not at promoters, and because the consensus is a union
  they propagate: 25,871 of 80,887 final clusters on that panel rested on those four
  libraries alone. The pipeline now reads the threshold back and gates on it, and reports
  it as a new `Log2FoldThreshold` column in both QC tables.
- `AlignedReads`, `BelowMapqReads` and `FilteredOutReads` in the stats and QC tables,
  between `TrimmedReads` and `TotalReads`. That gap is the largest single loss anywhere in
  the pipeline (a median 74.9% per library, 2,149 M of 2,954 M reads, on one 34-library
  panel) and nothing recorded where it went. Three things could account for it with
  opposite implications, and the QC could not tell them apart: reads that never aligned,
  reads below `filtering/alignment_mapq`, and reads dropped by `exclude_flags`,
  `min_alignment_length` or `max_mismatch`. The counts come from a new pass-through awk in
  the alignment pipe, since the unfiltered alignment is never written to disk and the
  flagstat the pipeline keeps runs on the filtered BAM, so it reports 100% mapped for
  every library. Full per-library detail, including a MAPQ histogram, lands in
  `qc/<sample>.aln.raw.txt`.
- `program/keep_unmapped_sample`, an integer, default 0: when above zero, keep that many
  unmapped reads per library as `qc/<sample>.unmapped.sample.fastq.gz`. With `keep_bam`
  and `keep_trimmed_fastq` off, nothing about the unaligned fraction survives a run, so a
  contaminant screen afterwards means aligning again.
- STAR's `Log.final.out` is kept rather than deleted, appended to
  `qc/<sample>.aligner.log`. STAR keeps unmapped reads out of its BAM entirely, so its own
  report is the only place its read accounting exists, and the new summary reads the input
  read count from it.
- A README caveat that every per-class QC figure (`PretRNAPct`, `miRNAPct`, any rRNA
  number) is computed from the filtered BAM and therefore measures the unique-mapping tail
  of those species rather than their abundance.
- `qc/min_cs_frip_final` and `qc/min_pct_nuclear_final`, gates for the final QC table.
  Both tables report `csFRiP`, but they measure it against different regions (the initial
  merged TSS set against `tss.final.bed`), so one number cannot serve both: every
  library's final `csFRiP` falls when the consensus shrinks, 0.030 to 0.058 on one
  22-library panel, so a gate tuned so the initial table flags exactly the four collapsed
  libraries went on to flag 8 of 22 in the final table, four of them sound. Unset by
  default, falling back to `min_cs_frip` and `min_pct_nuclear`, so existing configs
  behave exactly as before.
- `StatusReason` in both csRNA QC tables, naming the gates that tripped. `Status` on its
  own said only Ok or FAIL, which is opaque once more than one gate is in play.
- `StatusBasis` in both csRNA QC tables, `initial` or `final-advisory`, recording that
  `collect_consensus_tss` consumes the initial `Status` and nothing consumes the final
  one.

### Changed
- **`qc/min_log2_fold` changes which libraries pass QC**, and with
  `filtering/exclude_failed_from_consensus` on it changes the consensus set, so an
  existing analysis will shift. Compare before adopting, as with the `-gtf` entry in
  0.8.0. On the panel above it removes the four collapsed libraries and the 25,871
  clusters that rested on them alone.
- `rule align` runs a counter in the alignment pipe, so its shell command has changed and
  Snakemake will want to realign existing runs. The alignment itself is untouched: the awk
  passes the stream through byte for byte, which the test suite asserts.
- `tss/<sample>.stats.txt` is now post-processed by `find_tss_initial`, and is a declared
  output of that rule rather than an untracked side effect. HOMER writes
  `Fraction of stable transcript TSS clusters: 0.00%` and, with no annotation,
  `Fraction Promoter-Distal TSS clusters: 100.00%` whether or not the input those need
  was supplied, so both read as measurements; each is rewritten to `na` when HOMER's own
  report shows it had nothing to compute it from (an empty true-positive set on the
  relevant side). A real measurement is never overwritten, and anything parsing that file
  should expect `na` in those two positions.
- `collect_consensus_tss.R` sends its warnings to stderr rather than stdout: the notice
  that a library failed QC and still contributes to the union, and the one about clusters
  wider than their chromosome. On stdout they landed among the progress narration in the
  Snakemake log and were easy to miss in a long run.

### Fixed
- Editing `qc/min_cs_frip` or `qc/min_pct_nuclear` and rerunning was a silent no-op. Both
  reached the QC scripts through `snakemake@config`, which no rerun trigger can see, so
  the stale `Status` stayed on disk with nothing in the log to say the new gate had been
  ignored. Both now arrive as rule params, which `--rerun-triggers params` (in
  Snakemake's default trigger set) acts on. `collect_consensus_tss.R` had the same bug
  for `filtering/tss_min_reps`, reading the config copy while the rule was already
  passing it as a param; no `snakemake@config` read is left in any script.
- A QC gate that cannot be evaluated (an `NA` `csFRiP` from a library with no nuclear
  reads) previously left `Status` as `NA`, which `collect_consensus_tss` treated as
  neither Ok nor FAIL. It now counts as a failure, named in `StatusReason`.
- The README pointed at `qc/qc_final_cs.txt` and `qc/qc_final_in.txt`; both have been at
  the output root since the output paths were reorganised.

## [0.9.0] - 2026-08-04

### Removed
- `program/homer/tss/rnaseq_tagdir` and `rnaseq_stability_only`, so HOMER's `-rna`
  is no longer wired in. `-gtf` stays, and it carries the whole benefit anyway:
  measured on the test genome, `-gtf` alone yields the same true and false positive
  sets (15 and 5) and the same chosen input threshold (0.415) as `-gtf` with `-rna`,
  because the RNA-seq side only adds the stable-transcript filter on top. Dropping it
  removes an option that required users to build a HOMER tag directory by hand with a
  splice-aware aligner, since this pipeline's trim and alignment settings cannot make
  one correctly, along with a strand-orientation trap that would have made the
  stability call meaningless rather than merely noisy if got wrong.
- The end-to-end step and wiring assertions now cover `-gtf` only.

### Fixed
- The README no longer cites `Skipping TSS assignment (can't find file for genome ...)`
  as evidence that the annotation is missing. That line comes from the `annotatePeaks`
  calls inside HOMER wanting a full genome directory and appears whether or not `-gtf`
  was given; the usable signal is the TP and FP counts, and the documentation now shows
  the before and after for both.

## [0.8.1] - 2026-08-04

### Fixed
- `program/homer/tss/rnaseq_tagdir` and `gtf` were accepted without being checked
  beyond existence, so a path that was not a HOMER tag directory, or not a GTF,
  passed start-up and then failed deep inside HOMER. Both are now validated up
  front: the GTF for nine tab-separated fields, the tag directory for a
  `tagInfo.txt`, with errors that say what the option actually wants.

### Added
- `program/homer/tss/rnaseq_stability_only` passes HOMER's `-noFilterRNA`, so
  RNA-seq can be used for the stable-transcript columns without filtering TSSs.
  Supplying `rnaseq_tagdir` does filter by default, which is HOMER's behaviour and
  was not documented: on the synthetic data it took clusters passing the RNA check
  from 29 down to 17 passing both. Default False, matching HOMER.
- End-to-end coverage for both optional annotation inputs, plus a GTF fixture for
  the synthetic genome. `tests/e2e/check_annotation_wiring.py` asserts on HOMER's
  own report rather than on the flags being present, since both options are easy to
  wire up so that they look fine and do nothing. It checks that HOMER read the GTF,
  built non-empty TP and FP sets from it (TP 15, FP 5 on the fixture), chose the
  input threshold from the data instead of falling back to `-defaultLog2Fold`, and
  applied a real RNA-seq threshold. Writing it caught that the previous
  "Skipping TSS assignment" line is not a usable signal: it comes from the internal
  `annotatePeaks` calls wanting a full HOMER genome directory and appears whether or
  not `-gtf` was given.

### Changed
- Documented that `rnaseq_tagdir` is a tag directory rather than FASTQ, and why the
  pipeline cannot build it: RNA-seq reads are spliced and longer, while the trim
  step caps reads at `max_read_length` and the aligners are configured for short
  unspliced reads, so routing mRNA-seq through them would lose junction-spanning
  reads and truncate the rest.

## [0.8.0] - 2026-08-04

### Added
- `program/homer/tss/gtf`: optional gene annotation passed to HOMER as `-gtf`.
  Without it HOMER cannot build the annotated-TSS and exon sets it uses as true and
  false positives to choose each library's enrichment threshold, so it falls back to
  `-defaultLog2Fold` and writes the promoter-distal and stable-transcript columns as
  placeholders (100.00% / 0.00%). Optional and empty by default, because the pipeline
  must keep working without an annotation, but it is the largest single quality lever
  in the configuration. Setting it changes TSS calls, so compare against an existing
  consensus before adopting it. `program/homer/tss/rnaseq_tagdir` likewise passes
  `-rna`, enabling HOMER's stable-transcript filter.
- `program/homer/tss/program`: choice of HOMER caller. `findcsRNATSS.pl` announces
  itself as a legacy placeholder for `findcsRNATSR.pl`; the two are the same program
  (81 diff lines, all renames and help text) and their tables are identical, verified
  on a 34,891-cluster library with zero differing rows. The pipeline now normalises
  TSR's `.tsr.txt` / `.alltsr.txt` names to the `.tss.txt` the rest of it expects. The
  default stays on the legacy name because TSR ships an unpatched `use lib` pointing at
  its author's install, so it only resolves `HomerConfig` when HOMER's bin is on
  `PERL5LIB`.
- `program/homer/tss/pseudo_count`, `default_log2_fold` and `local_fold` expose HOMER's
  `-pseudoCount`, `-defaultLog2Fold` and `-L`, all defaulting to HOMER's own values.
  `pseudo_count` is what lets a cluster with zero input coverage pass the enrichment
  test automatically. `-cpu` is now passed so HOMER uses the threads the rule reserves.
- `filtering/exclude_failed_from_consensus` (default False, so behaviour is unchanged):
  drop csRNA libraries whose initial QC `Status` is FAIL from the consensus union. The
  union is a union, so a contaminated library's spurious clusters otherwise enter the
  shared set and everything downstream quantifies them. A warning naming the libraries
  is printed either way.
- `run_info.txt`: per-run provenance. Pipeline version and git commit (flagged when the
  working tree is dirty), host, resolved config, sample-table checksum, and the version
  of every tool the run used.
- Golden-value regression test. `tests/data/golden/qc_final_*.txt` are compared
  value-by-value after the default e2e run, so a metric changing silently now fails the
  suite; the structural assertions never noticed such changes. Refresh deliberately with
  `tests/scripts/update_golden.sh <e2e outdir>`. Skipped for the per-aligner runs, where
  a different aligner could legitimately shift values.

### Changed
- `tss_size_composition` is split into a per-library rule and a join, so Snakemake runs
  the tag-directory pass for each library in parallel instead of looping over them in
  one job.
- `tss.consensus.sizes.txt` now carries `.top1` … `.topN` up to the new
  `filtering/tss_top_sizes_max` (default 5) instead of a single `.topn`, so
  `tss_top_sizes_n` can be retuned anywhere in that range without re-reading the tag
  directories. Keeping the full per-cluster length histogram would serve any n exactly
  but runs to roughly 2.8M rows per library on a plant-sized dataset, against a handful
  of columns here.
- Documented what the QC table does not mean: that `csFRiP` and friends measure
  different regions in the initial and final tables, that the contamination percentages
  are signal-relative rather than library-relative, that a high duplicate rate can mean
  a *better* library because concentrated reads duplicate more, and that `csEnrichment`
  could not separate 1.97 from 2.00 at two to three libraries per arm.

## [0.7.0] - 2026-08-04

### Fixed
- `gather_stats.sh` summed organelle tag counts in bash integer arithmetic. HOMER
  writes tag counts as floats, and any fractional total (multi-mapper weighting,
  `bam2td -keepAll`) comes back from awk in `%.6g` scientific notation, which bash
  rejects with "syntax error: invalid arithmetic operator", failing the rule. The
  sum and rounding now happen inside awk.
- `collect_consensus_tss.R` could place clusters outside their chromosome:
  widening to 150 bp and the overlap-resolution shifts were unclamped and no
  seqlengths were set, so a TSS near a chromosome end could be written to BED with
  a negative start or a coordinate past the chromosome. Such clusters are now slid
  back inside, preserving their width, and a cluster wider than its own chromosome
  is left alone and reported.
- `collect_consensus_tss.R` picked one of each symmetric overlap pair with
  `seq(1, length(hits), by = 2)`, which assumes the two members of a pair land
  adjacent in the `Hits` object. That holds for equal-width sorted ranges but is
  not guaranteed, so it now selects on `queryHits < subjectHits`.
- `tss_size_composition.R` built its per-(cluster, length) key in integer
  arithmetic, which overflows silently to `NA` for a large genome sequenced with
  long reads (roughly a million clusters and multi-kilobase reads). Now double.
- A failed `findcsRNATSS.pl` call is no longer silent. The rule still writes an
  empty TSS set so one bad library cannot stall a run, but it now says so on
  stderr and in the log, because the only other symptom was the affected library
  quietly vanishing from the consensus union and a confusing `bedtools` error
  downstream.

### Changed
- `qc_final_cs.txt` and `qc_final_in.txt` no longer carry `FinalRiP`, `FinalFRiP`
  and `FinalEnrichment`. They were literal copies of `csRiP`, `csFRiP` and
  `csEnrichment` in the same table. `FinalTSSDetected` stays, since it is a
  distinct metric. Note the naming subtlety the duplication was hiding, now
  documented in the README and in the script: in the final tables `csRiP` /
  `csFRiP` / `csEnrichment` are computed against `tss.final.bed`, whereas the same
  names in the initial tables are computed against the initial merged set, so the
  two files must not be compared column-by-column.
- `miRNAPct` now uses the same denominator as `PretRNAPct`, contaminant reads over
  contaminant plus reads-in-TSS, instead of contaminant over nuclear reads. The two
  contamination percentages are reported side by side and were on different bases,
  differing by roughly the csFRiP factor. Values shift up by around 15% relative;
  ratios between libraries are essentially unchanged.
- Corrected the documented rationale for the small-RNA size filter (README,
  `tss_size_composition.R`, `normalize_tss_quantification.R`). It previously said
  enrichment over the input fails to reject siRNA clusters where the input is too
  shallow to measure a background. Measured properly, shallow input is a minor
  sub-mode: across the whole population of siRNA-dominated clusters only 1-2% have
  zero input coverage, and a 2-fold enrichment filter rejects 99% of them in a good
  library. What the filter really guards against is loss of library quality, since a
  library that has lost its capped signal is proportionally richer in siRNA than its
  own input: siRNA-dominated clusters clearing the enrichment filter run 1% in a good
  Arabidopsis library, 7% without the AP phosphatase and 61% in the worst library of
  the panel. Read size does not degrade that way. Documentation only, no behaviour
  change.

### Added
- Top-lengths filter, a second read-size composition filter that needs no prior
  knowledge of which lengths are contaminating. `tss_max_top_sizes_fraction`
  (default 1, disabled) drops a TSS cluster when more than that fraction of its
  reads fall in its `tss_top_sizes_n` (default 2) commonest read lengths.
  Genuine initiation is heterogeneous, spreading over tens of lengths, while a
  discretely processed RNA with fixed 5' and 3' ends puts nearly everything into
  one or two. Gated by the existing `tss_srna_min_reads` and
  `tss_srna_min_samples`, whose descriptions now say they govern both size
  filters.

  Measured on Arabidopsis csRNA-seq over 113,557 consensus clusters, `n = 2` with
  a 0.8 ceiling flags 0.09% of protein-coding clusters and 17% of transposon
  clusters. It catches contaminants that fall outside any plant small RNA size
  class, for example an abundant 27-28 nt 5'-polyphosphate species over a
  BRODYAGA1A element (top-2 share 0.83) that `tss_srna_sizes: [21..25]` misses
  entirely.
- `tss.consensus.sizes.txt` gains a `.topn` column per library. Columns are now
  written on demand: `.reads` always, `.srna` only when `tss_srna_sizes` is set,
  `.topn` only when the top-lengths filter is on. `tss_top_sizes_n` is baked into
  the table, so changing it re-runs `tss_size_composition`.
- Unit coverage for the new filter: two more `normalize` invocations (7 of 10
  clusters kept at `tss_srna_min_samples` 1, 8 at 2) and a third
  `size_composition` invocation checking the `.topn` column layout. The tag
  fixtures now give most clusters a spread of 5-10 read lengths so both filters
  have real negatives as well as positives.

### Changed
- `tss_srna_min_reads` default raised from 30 to 100. Size composition is
  estimated from few reads at weakly expressed clusters, so the low floor cost
  specificity: measured on Arabidopsis data the top-lengths flag rate runs 2.4%
  at 30-50 reads per cluster against 0.1-0.3% above 300. Genuine single-locus
  contaminants carry thousands of reads, so the higher floor keeps them. Only
  affects runs that enable a size filter; `tests/data/config.yaml` pins 30 so the
  unit fixtures still exercise the floor.

## [0.5.0] - 2026-08-04

### Added
- Read-size composition filter on the consensus TSS set, driven by four new
  `filtering` config keys: `tss_srna_sizes` (array of read lengths, e.g.
  `[21, 22, 23, 24, 25]`), `tss_max_srna_fraction` (default 0.5),
  `tss_srna_min_reads` (default 30) and `tss_srna_min_samples` (default 1).
  A TSS is dropped from `tss.final.bed` when more than `tss_max_srna_fraction`
  of its reads fall in the listed sizes, in at least `tss_srna_min_samples`
  csRNA libraries that have at least `tss_srna_min_reads` reads there.
  Defaults to `[]`, which disables the filter.

  This catches uncapped small RNAs that the existing filters cannot. In plants
  the 21-25 nt siRNAs, above all the 24 nt Pol IV class over transposons,
  survive TEX/AP well enough to be called as TSS clusters, and enrichment over
  the input cannot reject them wherever the input library is too shallow to
  measure a local background. Read length separates them from genuine
  initiation, which is not confined to one small size class.
- New rule `tss_size_composition` and script
  `workflow/scripts/tss_size_composition.R`, measuring per-cluster read counts
  and small-RNA-sized read counts from the HOMER tag directories into
  `tss.consensus.sizes.txt`. The length column of the tag files supplies the
  sizes, so no extra HOMER run is needed. When `tss_srna_sizes` is empty the
  rule writes the cluster list and does not read the tag directories at all.
- Unit test `tests/unit/size_composition.smk` with tag-directory fixtures, plus
  two further `normalize` invocations covering the new filter at
  `tss_srna_min_samples` 1 and 2.

### Changed
- `normalize_tss_quantification` now stops with an explanatory message when the
  filters leave no TSSs at all, rather than failing inside edgeR with
  "'counts' must contain at least one value".

### Fixed
- `tests/scripts/generate_unit_fixtures.py` writes the `RawReads` and
  `TrimmedReads` stats columns added in 0.3.0, so the `qc_initial` and
  `qc_final` unit tests pass again.
- `tests/run_tests.sh` propagates unit-test failures out of `run_unit_smk`
  instead of returning the status of its last conditional, which had been
  reporting failed assertions as passes.

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
