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

For a sound library the threshold HOMER chooses is a little stricter than the constant
it replaces, and TSS counts fall slightly. For a library that has lost its capped signal
it can come out at or below zero, and a threshold at or below zero accepts TSSs that
carry no more signal than their own input, and such a TSS cannot be genuine. Such a
library then reports many more TSSs with an annotation than without it, so the
annotation makes it look more productive rather than less.

The choice happens inside HOMER and there is no flag to put a floor under it, so the
pipeline reads the threshold back out of `tss/<sample>.stats.txt` and gates on it here.
Falling back to `default_log2_fold` when the key is unset is the defensible bound, since
a threshold chosen from the data should not end up looser than the constant it replaced.

## Why one failed library matters to everything

The consensus set is a union, so a library that fails QC does not simply carry poor
numbers of its own. The spurious TSSs it calls join the shared set, every other library
is then quantified over them, and a single bad library can contribute a large block of
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

Abundant uncapped small RNAs are the one contaminant the other filters cannot reach,
since they survive the enzymatic depletion well enough that HOMER calls them as TSSs.
Enrichment over the matched input holds them off while the csRNA library is clean, but
it weakens as the library degrades, because a library that has lost its capped signal is
proportionally richer in small RNA than its own input is. Read length does not degrade
that way, and genuine initiation is never confined to a single small size class, which
is what makes length the more durable signal of the two.

## Why there is a second, list-free size filter

`tss_srna_sizes` catches only contaminants whose length you can name in advance.
`tss_max_top_sizes_fraction` needs no list and asks instead how concentrated a TSS is in
its own commonest lengths. Genuine initiation is heterogeneous, because a promoter fires
across a window and the resulting RNAs vary in length, so a genuine TSS spreads over
tens of lengths. A discretely processed RNA has one 5′ end and one 3′ end and puts
nearly everything into one or two, which is how this filter catches species that fall
outside any small RNA size class.

A handful of reads gives a poor estimate of concentration, so the flag rate rises at
weakly expressed TSSs. That is why `tss_srna_min_reads` defaults to 100: a genuine
single-locus contaminant is usually very abundant, so a read floor keeps it while
dropping most of the small-sample noise. Lower the floor if you would rather be
sensitive.

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
positions alone would treat much of the genuine diversity as duplication. Two reads that
share both ends are far more likely to be copies of one molecule, which is what a
duplicate should mean.

Restricting the count to the final TSS set decides the answer rather than refining it. A
library made of scattered degradation background has more distinct positions genome-wide
than a good library does, so measured genome-wide the worst libraries can rank highest
on complexity. Restricted to called TSSs the ranking inverts to something usable.

The table reports Chao1 in place of a depth-matched count. To match depth across a set
of libraries you must cut every library down to the shallowest one, which either
discards most of a good library's data or refuses to answer for the libraries most worth
judging. Chao1 needs only singleton and doubleton counts, and it has a value for every
library. Good's coverage answers the separate question of whether more sequencing would
pay. Nothing gates on any of the three, partly because Chao1 still inflates for
background-heavy libraries and partly because saturation restates purity to a degree.

## Why `snRNA5pPct` is the one to watch and `tRNA5pPct` is not

Sm-class snRNAs are capped Pol II transcripts, so a working library puts nearly all of
their signal on the annotated 5′ end, and losing that precision is a direct symptom of a
phosphatase step that has failed. tRNA 5′ ends carry a monophosphate left by RNase P and
become ligation-competent when the same step fails, so tRNA precision moves the opposite
way, but only in the most extreme cases, which leaves it too insensitive to gate on. The
table carries it as a cross-check on which chemistry failed.

Do not combine the two into a ratio or a difference. The tRNA behaviour swamps the snRNA
signal, so either combination loses the clean separation the snRNA measure has on its
own.

## `findcsRNATSS.pl` against `findcsRNATSR.pl`

HOMER's `findcsRNATSS.pl` calls itself a legacy placeholder and points at
`findcsRNATSR.pl`, but the two are the same program: the diff is renames and help text,
and a row-for-row check on a full-size library found their output tables identical. The
default stays on the legacy name only because `findcsRNATSR.pl` ships with an unpatched
`use lib` pointing at its author's own install.

## Two metrics that do less than they appear to

`csEnrichment` does not discriminate at the group sizes these experiments usually have.
With two or three libraries per condition it has failed to separate conditions that
demonstrably differed, so power comparisons belong on the contamination and survival
metrics instead.

A duplicate rate is not a quality measure either. Duplication at matched depth is
complexity convolved with concentration, since reads pile up where the signal is, so a
library that puts more of its reads into TSSs duplicates more at the same depth. The
least duplicated library in a set can easily be the worst one, because its reads scatter
over degradation background. The complexity columns exist to answer this question
properly.
