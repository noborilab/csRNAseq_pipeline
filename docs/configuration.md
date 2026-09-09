# Configuration

Snakemake validates `config/config.yaml` against `workflow/schemas/config.schema.yaml`,
which carries a description and a default for every key. The entries below are the ones
where that description alone does not really tell you what to set.
[why-these-defaults.md](why-these-defaults.md) explains the reasoning behind the
defaults and the QC gates.

## `sample_table`

A TSV with at least three columns, `sample_name`, `sample_type` and `read_r1`. The
optional columns are `replicate`, `read_r2` and `input_name`.

- `sample_name`: base name for the library. Matched csRNA and input libraries share it
  by default, and so do replicates of the same condition. The pipeline groups libraries
  by this name during normalization.
- `sample_type`: `csrna` or `input`.
- `read_r1`: path or paths to R1 FASTQ files, comma-separated for several files
  belonging to one sample.
- `replicate`: integer replicate number, defaulting to 1 when the column is absent.
  Matched csRNA and input libraries must carry the same value.
- `read_r2`: R2 paths for paired-end data, left blank for single-end. A dataset may mix
  the two by leaving the column empty on the single-end rows.
- `input_name`: on csRNA rows only, the `sample_name` of the input library to pair with.
  Use it when several csRNA conditions share one input library; blank means the row
  pairs with its own `sample_name`.

```
sample_name	sample_type	read_r1	read_r2	input_name
sample1	csrna	sample1_1.csrna.r1.fq.gz,sample1_2.csrna.r1.fq.gz		
sample1	input	sample1.input.r1.fq.gz	sample1.input.r2.fq.gz	
sample2	csrna	sample2.csrna.r1.fq.gz		sample1
```

## `chrom_sizes`

Two columns, chromosome name and size in bp, with no header. A samtools fai index works
directly. This file decides which chromosomes survive the whole pipeline: only those
listed reach the merged TSSs, the consensus set and the bigWigs. Omit organellar and
other non-nuclear chromosomes here to exclude them.

## `qc / organelle_chroms`

Chromosomes to treat as organellar when computing `PctNuclear`. The pipeline subtracts
reads that map to them from the total aligned count. Default `["Pt", "Mt"]`, for
Arabidopsis.

```yaml
qc:
  organelle_chroms: ["chrM"]          # human / mouse
  organelle_chroms: []                # skip organelle subtraction
```

## `program / genome_index`

Path to the aligner's index. When it does not exist and you have set `genome_fasta`, the
pipeline builds it at `genome_index`. Use a directory for STAR and a file prefix for
the other aligners. Keep references immutable within a project; when changing a FASTA,
use a new index path and output directory.

## `program / homer / genome`

A HOMER genome name such as `hg38`, `mm10` or `tair10`, matching a genome configured in
your HOMER installation, or else a path to a FASTA of your own genome.

## `program / homer / tss / gtf`

Gene annotation passed to HOMER as `-gtf`. It is empty by default, and you only need to
set it when `program / homer / genome` points at a plain genome FASTA, because an
installed HOMER genome already carries an annotation of its own.

If HOMER has no annotation from either source, it has no annotated-TSS and exon sets to
use as true and false positives, so it cannot choose an enrichment threshold from the
data and falls back to `default_log2_fold` instead. It does so quietly: the log reports
empty TP and FP sets, and `tss/<sample>.stats.txt` reports the promoter-distal and
stable-transcript fractions as placeholders, which the pipeline rewrites to `na`. With a
GTF the same library reports non-empty sets and a threshold taken from its own data.

Supplying an annotation **changes your TSS calls**, so it is worth comparing a run that
has one against your existing consensus instead of switching over mid-project.

```yaml
program:
  homer:
    tss:
      gtf: "/path/to/annotation.gtf"
```

While judging whether the GTF took effect, ignore the `Skipping TSS assignment (can't
find file for genome ...)` line. It comes from the `annotatePeaks` calls inside HOMER
wanting a full HOMER genome directory for their own TSS-distance annotation, and it
appears either way. The TP and FP counts are the signal.

## `program / homer / tss / program`, `pseudo_count`, `default_log2_fold`, `local_fold`

`findcsRNATSS.pl` calls itself a legacy placeholder and points at `findcsRNATSR.pl`. The
two are the same program, and their output tables are identical row for row, so
switching changes only file names and labels, which the pipeline normalises. The default
stays on the legacy name because `findcsRNATSR.pl` ships with an unpatched `use lib`
pointing at its author's own install, so it finds `HomerConfig` only when HOMER's `bin`
is on `PERL5LIB`. Check that inside your container before setting `program:
"findcsRNATSR.pl"`.

Three HOMER thresholds are configurable, each defaulting to HOMER's own value:

- `pseudo_count` (1.0). HOMER adds this to both sides of the enrichment ratio, which is
  stabilizes ratios at low coverage. Increasing it reduces ratios driven by small
  counts, but neither its default nor a larger value requires nonzero input coverage
  or guarantees that a cluster passes the other calling criteria.
- `default_log2_fold` (1) is the threshold HOMER uses whenever it has no annotation.
- `local_fold` (2) is the local enrichment required during initial TSS identification.

## `filtering / exclude_failed_from_consensus`

The consensus set is a union, so a library that fails QC does not merely have poor
numbers of its own: its calls enter the shared set and everything downstream
quantifies them, including any false calls. A failed QC gate does not itself label
individual loci as false. Set this to drop failed libraries from the union. It
defaults to `False`, and the pipeline prints a warning that names the offending
libraries either way.

```yaml
filtering:
  exclude_failed_from_consensus: true
```

`Status` is FAIL when `csFRiP` is below `qc / min_cs_frip`, nuclear percentage is below
`qc / min_pct_nuclear`, or the enrichment threshold HOMER chose is below `qc /
min_log2_fold`, and `StatusReason` says which. All three gates reach `qc_initial_tss` as
rule params, so editing one and rerunning recomputes the table; Snakemake's default
`--rerun-triggers` set includes `params`. Consider `tss_min_reps: 2` alongside it, so no
single library can inject a TSS on its own.

## `qc / min_log2_fold`

The lowest enrichment threshold a library's TSS calling may use before the library fails
QC. Unset by default, in which case it falls back to `program / homer / tss /
default_log2_fold`.

With an annotation HOMER chooses this threshold per library from annotated TSS and
exon distributions. A low value can flag weak separation between initiation signal and
background, but also depends on annotation quality, the input library, and the biology.
It is not proof that a library failed or that every accepted locus is false. HOMER offers
no flag to floor the selected threshold, so this pipeline can gate the library here.

The fallback to `default_log2_fold` is a conservative policy choice, not a statistical
requirement. Inspect the paired libraries and annotation before excluding a sample.
Set this key explicitly to use a separately calibrated bound.

```yaml
qc:
  min_log2_fold: 1.0
```

A library that trips this gate is FAIL with `Log2FoldThreshold` in its `StatusReason`,
and `filtering / exclude_failed_from_consensus` is what keeps its TSSs out of the union.

## `program / keep_unmapped_sample` and `program / keep_below_mapq_sample`

How many discarded reads to keep per library, from each of the two classes that account
for most of the alignment loss. Both default to 0, which writes nothing.

| Key | File | What it holds |
|---|---|---|
| `keep_unmapped_sample` | `qc/{sample}.unmapped.sample.fastq.gz` | reads that never aligned |
| `keep_below_mapq_sample` | `qc/{sample}.belowmapq.sample.fastq.gz` | reads that aligned but fell below `filtering / alignment_mapq` |

A large share of a csRNA-seq library can be lost at alignment, and with `keep_bam`,
`keep_trimmed_fastq` and `bwa_aln / keep_sai` all false neither the aligned nor the
unaligned reads survive the run, so a contaminant screen afterwards means aligning
again. A few hundred thousand reads of each costs a few MB.

```yaml
program:
  keep_unmapped_sample: 200000
  keep_below_mapq_sample: 200000
```

They are separate keys and separate files because the two classes mean opposite things.
A read that never aligned came from something that is not in your genome, or is too
short or too poor to place. A read below the MAPQ floor placed fine and placed somewhere
else equally well, which with bwa means MAPQ 0 and a multi-mapper. `BelowMapqReads` is
usually the larger class and the one to look at first when asking where a big loss went.

There are two further points about these files. Reads in the below-MAPQ sample are
written as sequenced, so a reverse-strand alignment is reverse-complemented back out of
alignment orientation first. And each sample is the **first** N reads of its class in
the order the aligner emitted them, which is roughly flowcell order rather than a random
draw, so it will tell you what these reads are but is no basis for quantifying anything.
A library with no reads in a class gets no file, and `qc/{sample}.aln.raw.txt` records
the zero.

With STAR the below-MAPQ sample comes from the alignment stream as with every other
aligner, while the unmapped sample needs `--outReadsUnmapped Fastx`, which the pipeline
adds.

## `qc / min_cs_frip_final` and `qc / min_pct_nuclear_final`

The initial and final QC tables both report `csFRiP`, measured against the initial
merged TSS set in one and against `tss.final.bed` in the other, so one number cannot
gate both. Every library's `csFRiP` falls when the consensus set shrinks, purely from
the set getting smaller, so a gate tuned on the initial table flags sound libraries in
the final one.

Both keys are unset by default and fall back to `min_cs_frip` and `min_pct_nuclear`.

```yaml
qc:
  min_cs_frip: 0.80         # gates the initial table, and the consensus through it
  min_cs_frip_final: 0.70   # gates the final table only
```

Only the initial `Status` feeds `collect_consensus_tss`. See the closing note in
[qc-metrics.md](qc-metrics.md).

## `filtering / tss_min_cpm` and `filtering / tss_min_samples`

The pipeline applies this filter after quantification and before TMM normalization. A
TSS reaches `tss.final.bed` only if its CPM, counts per million library-size reads
without TMM, is at least `tss_min_cpm` in at least `tss_min_samples` csRNA samples. The
defaults of 0 and 1 keep every TSS. Here the library size is the sum of counts
in all consensus TSSs before filtering, not all sequenced or mapped reads.

```yaml
filtering:
  tss_min_cpm: 1
  tss_min_samples: 2
```

## `filtering / tss_srna_sizes` and the read-size composition filter

Abundant processed small RNAs can survive library preparation and be called as TSSs.
Concentration at 21–25 nt can flag plant small-RNA contamination, but read length alone
does not establish RNA origin or cap chemistry. Size selection, trimming, sequencing
length and sampling depth also shape this distribution. Validate the filter against
known initiation loci and examine the clusters it removes.

List the contaminating lengths in `tss_srna_sizes`. The pipeline then drops a TSS when
more than `tss_max_srna_fraction` of its reads fall in those sizes. Only libraries with
at least `tss_srna_min_reads` reads in the TSS get a vote, and `tss_srna_min_samples` of
them must agree before the TSS goes. The default `[]` disables the filter, and while the
filter is off, `tss_size_composition` does not read the tag directories at all.

```yaml
filtering:
  tss_srna_sizes: [21, 22, 23, 24, 25]
  tss_max_srna_fraction: 0.5
  tss_srna_min_reads: 100
  tss_srna_min_samples: 1
```

`tss_srna_min_samples: 1` drops a TSS that any single library flags, which is the
sensitive choice when libraries differ in size selection: a 40-80 nt library cannot see
a 24 nt population that a 20-80 nt library shows plainly, so demanding agreement would
let the contaminant through. Raise it to keep a TSS that only one library objects to.

The per-TSS counts land in `tss.consensus.sizes.txt` whether or not the filter removes
anything.

## `filtering / tss_max_top_sizes_fraction` and the top-lengths filter

The size list above catches only contaminants whose length you can name in advance. This
filter needs no list and asks how concentrated a cluster is in its commonest read
lengths. Processed RNAs can be concentrated in one or two lengths, whereas many
promoter-derived libraries have broader distributions. This is a heuristic: a narrow
insert distribution is not proof of contamination, and a focused 5′ initiation site
can still produce heterogeneous RNA lengths. Technical size selection and truncation
must be considered before applying the filter.

Set `tss_max_top_sizes_fraction` below 1 to enable it and `tss_top_sizes_n` for how many
lengths to add up. The `tss_srna_min_reads` and `tss_srna_min_samples` gates apply here
too.

```yaml
filtering:
  tss_top_sizes_n: 2
  tss_max_top_sizes_fraction: 0.8
```

At a weakly expressed TSS the pipeline estimates concentration from very few reads, so
the flag rate depends on depth. That is why `tss_srna_min_reads` defaults to 100. Real
single-locus contaminants are usually very abundant, so the floor keeps them while
dropping most of the small-sample noise. Lower it to around 30 to be more sensitive.

`tss_size_composition` precomputes the share held by the 1 to `tss_top_sizes_max`
commonest lengths, 5 by default, and `normalize_tss_quantification` picks out the one
`tss_top_sizes_n` asks for, so retuning `tss_top_sizes_n` within that bound costs
nothing and it must not exceed it. Only `tss_srna_sizes` and `tss_top_sizes_max` force a
re-read of the tag directories.

## Stricter discovery options

The defaults retain a permissive discovery catalogue. For a study with at least two
biological replicates per condition, the following is a starting point to evaluate,
not a universally validated preset. Merge these keys into your existing configuration:

```yaml
filtering:
  exclude_failed_from_consensus: true
  tss_min_reps: 2
  tss_min_cpm: 1
  tss_min_samples: 2
  # Optional plant-specific screens; inspect losses at known promoters first.
  tss_srna_sizes: [21, 22, 23, 24, 25]
  tss_max_srna_fraction: 0.5
  tss_top_sizes_n: 2
  tss_max_top_sizes_fraction: 0.8
  tss_srna_min_reads: 100
  tss_srna_min_samples: 1
qc:
  min_cs_frip: 0.9
  min_pct_nuclear: 90
  min_log2_fold: 1
```

Calibrate the three QC gates against controls before using them to exclude libraries;
these are the existing default gates, not externally validated quality cutoffs.
`tss_min_reps` requires support within a condition, after QC exclusion. A condition
with fewer surviving libraries contributes no new loci, though its libraries are still
quantified at loci contributed by other conditions. QC exclusion affects discovery
only: failed libraries remain in count matrices, normalization and size-filter voting.
For a final biological analysis, remove rejected libraries from the sample table and
rerun in a new output directory; otherwise they can still influence those steps.

The CPM sample count applies across all csRNA libraries. The size filters can delete a
cluster globally because of one library (`tss_srna_min_samples: 1`); raising that to 2
reduces this risk but can miss contaminants only visible in one size-selection arm.
Low-count clusters below `tss_srna_min_reads` escape these size screens. Compare the
filters separately before combining them, report their losses, and avoid changing
thresholds to maximize a desired treatment difference. Final QC gates are descriptive
and may need separate calibration because their TSS regions differ.

## Optional csRNA/input ratio screen

`filtering / tss_min_cs_in_ratio` defaults to 0 (disabled). When enabled, the pipeline
quantifies input libraries over the exact same consensus regions as csRNA libraries,
including their adjusted boundaries. It retains a cluster if at least
`tss_min_ratio_samples` paired comparisons meet the requested ratio:

```
(csRNA_count / csRNA_non_organelle_depth) /
(max(input_count, 1) / input_non_organelle_depth)
```

Depth is `TotalReads - OrganelleReads` from the filtered tag directories, with organelles
defined by `qc / organelle_chroms`. These depths include non-organelle contigs omitted
from `chrom_sizes`, so use consistent references and organelle definitions. The input
floor is one read before depth normalization. This is an enrichment screen, not a
significance test; low counts remain uncertain. A threshold of 2 means twofold
normalized enrichment. It is optional even in the stricter configuration because HOMER
already screens enrichment, and this extra screen can remove weak genuine initiation.

Versions through 0.16.2 used raw count ratios and joined input counts from the initial
region set by reassigned IDs. Results with this option enabled must be regenerated;
old numerical thresholds do not have the same meaning.

## Annotation masks and biological scope

Providing `qc / mirnas` or `qc / trnas` also removes consensus regions overlapping
those annotations on either strand. This can remove primary miRNA initiation or nearby
and antisense promoters as well as processed-RNA contamination. Leave a mask empty if
those loci are in scope, and evaluate contamination using other QC evidence. The current
configuration couples these QC annotations to masking; it does not distinguish mature
processed RNA from primary transcription simply by overlap. The [original csRNA-seq
study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6836739/) includes primary miRNA initiation.
