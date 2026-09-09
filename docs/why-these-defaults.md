# Why these defaults

This file explains why the defaults and the QC gates have the values they do, which is
the one thing neither the schema nor the config comments have room for. It deliberately
carries no measurements, so nothing here depends on a particular dataset. The release
entry in `CHANGELOG.md` records what prompted each choice at the time.

## Why an annotation matters so much

Without an annotation, whether from an installed HOMER genome or from `program / homer /
tss / gtf`, HOMER has no annotated-TSS set and no exon set to use as true and false
positives, so it cannot choose an enrichment threshold from the data. It falls back to
`default_log2_fold` and says so only in passing:

```
Total TP (tss) regions: 0     Total FP (exon) regions: 0
Maximum CDF difference: -10000000000
Using the default input threshold (1)
```

With an annotation the same library reports non-empty sets and a threshold taken from
its own data. Running the pipeline's synthetic test genome, for instance, gives:

```
Custom annotation GTF file: .../annotation.gtf (using transcript_id)
Total TP (tss) regions: 15     Total FP (exon) regions: 5
Maximum CDF difference: 0.8
log2 fold vs. input: 0.415031488137169
```

While checking whether an annotation took effect, ignore the `Skipping TSS assignment
(can't find file for genome ...)` line, which comes from the `annotatePeaks` calls
inside HOMER wanting a full HOMER genome directory for their own TSS-distance annotation
and appears either way. The TP and FP counts are the signal.

## Why `qc / min_log2_fold` exists

A low annotation-derived threshold can flag weak separation of initiation signal from
background. It can also reflect the annotation, input composition or biological context;
it does not prove that every accepted locus is false. The threshold belongs to the
library-level calling procedure and is not the enrichment of every individual cluster.

The pipeline reports this threshold and can fail a library under a configured policy.
Using `default_log2_fold` as the fallback floor is conservative, not a statistical rule
that learned thresholds must be stricter than fixed ones. Calibrate against controls and
inspect other QC before excluding a library.

## Why one failed library matters to everything

The consensus set is a union, so a library that fails QC does not simply carry poor
numbers of its own. Its calls join the shared set, every other library is then
quantified over them, including any false calls, and a single bad library can contribute a large block of
the final set on its own. That is what `exclude_failed_from_consensus` is for, and
`tss_min_reps: 2` is worth considering alongside it so that no one library can inject a
TSS that no other library supports.

## Why the final QC gates are separate from the initial ones

Both QC tables report `csFRiP`, but the initial table measures it against the merged
initial TSS set and the final table against `tss.final.bed`. Every library's `csFRiP`
falls when the set behind it shrinks, purely because that set is smaller, so a threshold
tuned to flag exactly the collapsed libraries in the initial table will flag sound
libraries in the final one. `min_cs_frip_final` and `min_pct_nuclear_final` exist so you
can set the two independently, and each falls back to its initial counterpart when you
leave it unset.

## Why keeping samples of discarded reads is worth the disk

A large share of a csRNA-seq library can be lost at alignment, and with `keep_bam`,
`keep_trimmed_fastq` and `bwa_aln / keep_sai` all false, neither the aligned nor the
unaligned reads survive the run. Answering "what were those reads" afterwards then means
aligning everything again.

The pipeline keeps the two classes apart because they mean opposite things. A read that
never aligned came from something absent from your genome, or is too short or too poor
to place. A read below the MAPQ floor placed perfectly well and placed somewhere else
equally well, which with bwa means MAPQ 0 and a multi-mapper. The below-MAPQ class is
usually the larger of the two, and repetitive families are the usual explanation:
assemblies commonly represent the rDNA arrays as many near-identical copies, so an
aligner cannot place rRNA reads uniquely, and they never reach the filtered BAM. That
also explains why an rRNA percentage computed from the filtered BAM can come out
implausibly low.

## Why the size filter uses read length

Processed small RNAs can survive the preparation and contribute apparent initiation
clusters. Concentration in expected small-RNA lengths is a useful screen, especially
when csRNA/input enrichment becomes less informative in contaminated libraries. It is
not a direct assay of cap state or RNA origin. Size selection and trimming can narrow
the observed lengths of genuine initiation products too.

## Why there is a second, list-free size filter

The top-lengths filter can flag concentrated RNA populations outside a specified size
list. Many promoter-derived libraries span a range of lengths, but this is an empirical
pattern rather than a universal property of initiation. A focused start site and a
narrow insert-length distribution are different measurements.

The default floor of 100 reads limits unstable estimates at low depth. It also leaves
weak contaminants unfiltered. Both size filters remain disabled by default; evaluate
sensitivity and losses at independently supported initiation sites before enabling them.

## What the read-length columns can and cannot show

`MedianReadLength` is the stable summary. `P20ReadLength` and `P80ReadLength` describe
the shape of the distribution but do not separate library types, because size selection
moves every percentile together, and neither does the spread between them. A csRNA
library compared against its own matched input is the comparison that works.

`ModeReadLength` says little alone, since a healthy csRNA library is nearly flat in
length and the winning bin flips between replicates of the same condition.
`ModeReadFraction`, the height of that bin, is the one that separates, because a library
made largely of one processed species concentrates in a few lengths while initiation
spread over thousands of promoters does not. Nothing gates on it, both because it is a
shape statistic rather than a measurement of purity and because it largely restates
`PctNuclear` and `csFRiP`.

## Why the complexity columns count what they count

A fragment is a `(position, strand, length)` triple rather than a 5′ position alone,
because a single in-TSS position carries several distinct lengths, and a count of
positions alone would treat much of the genuine diversity as duplication. Identical coordinates can arise from independent molecules as well as PCR copies.
Without UMIs these are fragment-pattern diversity estimates, not molecule counts.

Restricting the count to the final TSS set decides the answer rather than refining it. A
library made of scattered degradation background has more distinct positions genome-wide
than a good library does, so measured genome-wide the worst libraries can rank highest
on complexity. Restricted to called TSSs the ranking inverts to something usable.

The table reports Chao1 in place of a depth-matched count. To match depth across a set
of libraries you must cut every library down to the shallowest one, which either
discards most of a good library's data or refuses to answer for the libraries most worth
judging. Chao1 uses singleton and doubleton fragment-pattern counts to estimate unseen
patterns under sampling assumptions. Good's coverage estimates the chance of seeing a
previously observed pattern, not the chance of discovering a new molecule or locus.
Neither replaces matched-depth comparisons or a library-complexity experiment. Nothing gates on any of the three, partly because Chao1 still inflates for
background-heavy libraries and partly because saturation restates purity to a degree.

## Why `snRNA5pPct` is the one to watch and `tRNA5pPct` is not

Sm-class snRNAs are capped Pol II transcripts, making annotated 5′-end concentration
a useful QC cross-check. Reduced precision is compatible with degradation background
or failed depletion, but can also reflect annotation and RNA processing. It is not a
specific biochemical diagnosis of phosphatase failure. tRNA-end signal can provide
complementary evidence, subject to mature versus precursor annotation and mapping bias.

Interpret both alongside matched input and preparation controls. Neither metric gates
the workflow, and their diagnostic value must be established in the organism and
library preparation being studied.

## `findcsRNATSS.pl` against `findcsRNATSR.pl`

HOMER's `findcsRNATSS.pl` calls itself a legacy placeholder and points at
`findcsRNATSR.pl`, but the two are the same program: the diff is renames and help text,
and a row-for-row check on a full-size library found their output tables identical. The
default stays on the legacy name only because `findcsRNATSR.pl` ships with an unpatched
`use lib` pointing at its author's own install.

## Two metrics that do less than they appear to

`csEnrichment` can have limited discriminatory power in small panels. A failure to
separate a particular set of conditions does not establish that the metric is generally
uninformative. Examine paired effects, uncertainty, contamination and alignment survival
together; none alone establishes assay specificity.

A duplicate rate is not a quality measure either. Duplication at matched depth is
complexity convolved with concentration, since reads pile up where the signal is, so a
library that puts more of its reads into TSSs duplicates more at the same depth. The
least duplicated library in a set can easily be the worst one, because its reads scatter
over degradation background. The complexity columns exist to answer this question
properly.
