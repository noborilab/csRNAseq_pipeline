# csRNA-seq pipeline

Author: Benjamin Jean-Marie Tremblay

This pipeline takes paired csRNA-seq (capped small RNA) and input (total RNA) libraries
and uses them to find and quantify transcription start sites (TSSs) across the genome. In
principle it should work with any organism for which you have a genome FASTA.

## Installation

Two environment specs live under `workflow/envs/`:

| File | Purpose |
|------|---------|
| `minimal_env.yaml` | Running the pipeline with bwa |
| `full_env.yaml` | The same plus every aligner, the test suite and development tools |

```bash
conda env create -f workflow/envs/minimal_env.yaml
conda activate csRNAseq
```

Three dependencies are not on conda, or not for every platform, and need installing by
hand:

- **bfqutils**: https://github.com/noborilab/bfqutils
- **bam2td**: https://github.com/noborilab/bam2td
- **HOMER**: on anything other than Linux, follow
  http://homer.ucsd.edu/homer/introduction/install.html, deleting the HOMER line from the
  yaml first if you still want to use the environment files.

Installing a HOMER genome with `perl configureHomer.pl -install tair10` is optional, and
convenient when one exists for your organism, since those genomes come with their own
annotation. If there is none, or you would rather not install it, point
`program / homer / genome` at a genome FASTA instead and supply the annotation yourself
through `program / homer / tss / gtf`.

## Running

With the environment set up and your config and sample table filled in, running the
pipeline is just an ordinary Snakemake invocation. You will want at least three cores,
since a couple of the HOMER programs always run three threads internally.

```bash
snakemake --cores 6 --configfile config/config.yaml        # add -n for a dry run first
```

## Testing

There is a test suite in `tests/`, built around a small synthetic genome and a set of
synthetic FASTQs. It covers schema validation, the correctness of the Snakemake DAG, the R
scripts on their own, and a full end-to-end run.

```bash
./tests/run_tests.sh all         # every layer
./tests/run_tests.sh schema      # schema validation only, and the only layer that
                                 # needs no biology tools installed
```

The layers can also be run individually as `validate`, `unit` and `e2e`. See
[`tests/README.md`](tests/README.md) for the full instructions, including how to
regenerate the fixtures.

## Configuration

`config/config.yaml` is validated against `workflow/schemas/config.schema.yaml`, which
carries a description and a default for every key. The entries below are the ones worth
understanding before a first run, and [`docs/configuration.md`](docs/configuration.md)
covers each of them properly.

| Key | What it does | Default |
|---|---|---|
| `sample_table` | TSV of libraries, keyed on `sample_name` and `sample_type`, with `input_name` to share one input library across conditions | required |
| `chrom_sizes` | Two-column TSV that decides which chromosomes survive the whole pipeline. Omit organellar chromosomes here to exclude them | required |
| `program / genome_index` | Aligner index, built from `genome_fasta` when it is absent | required |
| `program / homer / genome` | An installed HOMER genome name, or a path to your own genome FASTA | required |
| `program / homer / tss / gtf` | Gene annotation for HOMER, which you only need to supply when `program / homer / genome` is a plain FASTA, since installed HOMER genomes bring their own. With no annotation from either source, HOMER cannot pick a per-library enrichment threshold and falls back silently to `default_log2_fold`. Adding one changes your TSS calls, so it is better not to switch mid-project | empty |
| `program / homer / tss / pseudo_count`, `default_log2_fold`, `local_fold` | HOMER's own thresholds, exposed | HOMER's values |
| `qc / organelle_chroms` | Chromosomes treated as organellar when computing `PctNuclear` | `["Pt", "Mt"]` |
| `qc / min_cs_frip`, `min_pct_nuclear`, `min_log2_fold` | Gates behind `Status` in the initial QC table, which is the table that can change a run | 0.9, 90, unset |
| `qc / min_cs_frip_final`, `min_pct_nuclear_final` | Gates for the final QC table, which measures a different region and so needs its own | fall back to the above |
| `filtering / exclude_failed_from_consensus` | Keep a failed library's clusters out of the consensus union | `False` |
| `filtering / tss_min_reps` | Replicates a TSS must appear in to enter the consensus set | 1 |
| `filtering / tss_min_cpm`, `tss_min_samples` | CPM filter applied before TMM normalization | 0, 1 |
| `filtering / tss_srna_sizes` | Read lengths to treat as uncapped small RNA, dropping clusters dominated by them | `[]` |
| `filtering / tss_max_top_sizes_fraction` | Drops clusters concentrated in their own commonest read lengths, catching contaminants no size list names | 1, disabled |
| `program / keep_unmapped_sample`, `keep_below_mapq_sample` | Retain this many discarded reads per library, for asking afterwards what they were | 0 |

## Workflow

| Rule | Step |
|------|------|
| `build_index` | Build the aligner index from `genome_fasta`, when none is found at `genome_index`. |
| `trim` | Adapter trimming and length truncation, with bfqutils. |
| `align` | Alignment with bwa-aln, bwa-mem, STAR, bowtie2 or hisat2, then MAPQ and length filtering. Counts the unfiltered alignment on the way past into `qc/{sample}.aln.raw.txt`. |
| `make_tagdir` | Tag directory from the BAM, with bam2td. |
| `find_tss_initial` | Per-library TSS calling with HOMER, using the paired input library as background. Runs on input libraries too, whose calls are used for QC and small-RNA masking. |
| `make_raw_bedgraph` | Strand-separated raw bedGraphs, with HOMER `makeUCSCfile`. |
| `gather_stats` | Read counts, organellar counts, read-length summary and tag frequencies from the tag directories, plus the alignment loss from the per-library summaries. |
| `merge_initial_tss` | Strand-aware merge of the per-library TSS BEDs, dropping chromosomes absent from `chrom_sizes`. |
| `quantify_initial_{cs,in}_tss` | HOMER quantification of the merged sets in every library. |
| `qc_initial_tss` | The initial QC table, and the `Status` gates that can change the run. |
| `collect_consensus_tss` | The consensus set, requiring detection in `tss_min_reps` replicates. Narrow TSSs are padded, overlapping ones split, and those over miRNA or pre-tRNA loci removed. |
| `quantify_final_tss` | HOMER quantification of the consensus set in the csRNA libraries. |
| `tss_size_composition_library`, `tss_size_composition` | Per-cluster read-length composition, one job per library, joined into `tss.consensus.sizes.txt`. Reads nothing while both size filters are off. |
| `normalize_tss_quantification` | The CPM filter, then the two size-composition filters, then TMM normalization with edgeR. Writes `tss.final.bed`. Every filter defaults to a no-op. |
| `quantify_final_in_tss` | HOMER quantification of `tss.final.bed` in the input libraries. |
| `qc_five_prime`, `qc_complexity` | Per-library 5′-end precision and library complexity, feeding the final QC table. |
| `qc_final_tss` | The final QC table, re-derived from `tss.final.bed` rather than copied from the initial one. |
| `generate_normalized_bw` | Raw bedGraphs scaled to RPM and written as bigWig, with small RNA regions optionally masked. |
| `run_info` | Provenance: pipeline version and commit, host, resolved config, sample-table checksum, and the version of every tool the run used. |

## Output

Everything lands in `files / output_dir`, `results/` by default. The QC tables have a file
of their own, [`docs/qc-metrics.md`](docs/qc-metrics.md), which also covers the supporting
files under `qc/`.

| File | Contents |
|------|----------|
| `tss.final.bed` | The filtered consensus TSS set, BED6. Equal to `tss.consensus.bed` while the filters are at their defaults. |
| `tss.final.raw.txt` | Raw tag counts per TSS per csRNA library. Rows follow `tss.final.bed`. |
| `tss.final.cpm.txt` | TMM-normalized CPM counts, same row order. |
| `tss.final.in.raw.txt` | Raw counts for the same set in the input libraries, used by the final QC. |
| `norm_factors.txt` | edgeR TMM factors and RPM multipliers. |
| `tss.consensus.bed` | The consensus set before filtering. |
| `tss.consensus.sizes.txt` | Per library and cluster: total reads, reads at `tss_srna_sizes` lengths, and the share held by the 1 to `tss_top_sizes_max` commonest lengths. Holds cluster names alone while both size filters are off, since nothing then reads the tag directories. |
| `bw/{sample}.rpm.{pos,neg}.bw` | RPM-normalized bigWigs, one per strand. Reverse-strand scores are negative. |
| `qc_final_{cs,in}.txt` | Per-library QC on the filtered set. |
| `qc/qc_initial_{cs,in}.txt` | Per-library QC on the merged initial set, and the source of the `Status` that gates the consensus. |
| `run_info.txt` | Provenance for the run. |

## Reproducibility

Two runs on the same data write the same bytes. Everything HOMER produces is sorted on a
coordinate key on the way out of the rule that makes it, because `annotatePeaks.pl` and
`findcsRNATSS.pl` both emit rows in the order their threads finish, and `makeUCSCfile` is
given an explicit `-color` because it picks a random one otherwise.

Four things still differ between two runs, none of them a result: the `*.benchmark.txt`
wall-clock records, the tool logs, `run_info.txt`, which is meant to be unique per run,
and one line of `tss/{sample}.stats.txt`, where HOMER records the random draw it used for
its own temporary filenames. Run into a different output directory and a further set of
files differ only where a tool records its own command line, absolute paths included.

## Documentation

- [`docs/configuration.md`](docs/configuration.md), the config keys that need more than a
  schema description.
- [`docs/qc-metrics.md`](docs/qc-metrics.md), every QC column, what the numbers do not
  mean, and the supporting files under `qc/`.
- [`docs/why-these-defaults.md`](docs/why-these-defaults.md), the measurements behind the
  defaults and gates.
- [`tests/README.md`](tests/README.md), the test suite and its fixtures.
- `CHANGELOG.md`, what changed and why, release by release.

## Licence

GNU General Public License v3.0. See [`LICENSE`](LICENSE) for the full text.
