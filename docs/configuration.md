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
pipeline builds it under `intermediate_dir`. Do not point `genome_index` and
`genome_fasta` at the same path, because the pipeline then assumes the index already
exists.

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
  what lets a TSS with no input coverage pass automatically, since `log2((cs+1)/(0+1))`
  is large by construction. Raise it to demand real input coverage.
- `default_log2_fold` (1) is the threshold HOMER uses whenever it has no annotation.
- `local_fold` (2) is the local enrichment required during initial TSS identification.

## `filtering / exclude_failed_from_consensus`

The consensus set is a union, so a library that fails QC does not merely have poor
numbers of its own: the spurious TSSs it calls enter the shared set and everything
downstream quantifies them. Set this to drop failed libraries from the union. It
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

With an annotation HOMER chooses this threshold per library from the data and writes it
into `tss/<sample>.stats.txt` as `log2 fold vs. input`. Sound libraries get sensible
values. Libraries that have lost their capped signal get thresholds at or below zero,
which accept TSSs that carry less signal than their own input, and because the consensus
is a union those calls propagate to everything downstream. HOMER offers no flag to floor
the threshold, so the pipeline reads it back and gates on it here.

Falling back to `default_log2_fold` is the more defensible bound, because a threshold
chosen from the data should not end up looser than the constant it replaced. You can set
this key explicitly if you would rather gate independently of it.

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
defaults of 0 and 1 keep every TSS.

```yaml
filtering:
  tss_min_cpm: 1
  tss_min_samples: 2
```

## `filtering / tss_srna_sizes` and the read-size composition filter

Abundant uncapped small RNAs are the one contaminant the other filters cannot reach. In
plants the 21-25 nt siRNAs, above all the 24 nt Pol IV class over transposons, survive
the enzymatic depletion well enough that HOMER calls them as TSSs. Enrichment over the
input holds them off while the csRNA library is clean but degrades as the library
degrades, because a library that has lost its capped signal is proportionally richer in
siRNA than its own input. Read length is the signal that does not degrade, since genuine
initiation is never confined to one small size class.

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
filter needs no list and asks instead how concentrated a TSS is in its own commonest
read lengths. Genuine initiation is heterogeneous, because a promoter fires across a
window and the resulting RNAs vary in length, so a protein-coding TSS spreads over tens
of lengths. A discretely processed RNA has one 5′ end and one 3′ end and puts nearly
everything into one or two.

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
