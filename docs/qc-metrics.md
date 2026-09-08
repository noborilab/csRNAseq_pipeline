# QC metrics

`qc/qc_initial_cs.txt` and `qc/qc_initial_in.txt` carry the columns below, measured
against the initial merged TSS set. `qc_final_cs.txt` and `qc_final_in.txt` carry the
same columns plus `FinalTSSDetected`, `NConsensusTSS`, `NFinalTSS` and `NFilteredTSS`,
with every value re-derived from `tss.final.bed` rather than copied across. A column
such as `csFRiP` therefore measures a different region in each file, so the two are not
comparable row by row. Before acting on any of these numbers, it is worth reading
[Reading these numbers](#reading-these-numbers) below.

## Columns

| Column | Description |
|--------|-------------|
| `Log2FoldThreshold` | The enrichment threshold this library's TSS calling used, read back from `tss/<sample>.stats.txt`. Chosen by HOMER from the data when an annotation is available, otherwise `program / homer / tss / default_log2_fold`. Gated by `qc / min_log2_fold`. It reads `NA` when HOMER wrote no threshold at all. That happens for a library with no valid TSSs, where HOMER hits a division by zero on the way to the promoter-distal fraction and stops its report early. `NA` counts as a failure, and input libraries commonly show it. |
| `AlignedReads` | Reads with a primary alignment, counted before filtering. `TrimmedReads` minus this is the reads that never aligned. |
| `BelowMapqReads` | Aligned reads discarded by `filtering / alignment_mapq`, which in practice means multi-mappers. |
| `FilteredOutReads` | Aligned reads that cleared MAPQ and were then discarded by `exclude_flags`, `include_flags`, `min_alignment_length` or `max_mismatch`. |
| `MedianReadLength` | Median length of the reads that survived into the tag directory, from HOMER's own length histogram. Describes the filtered library, not the trimmer's output. |
| `P20ReadLength`, `P80ReadLength` | 20th and 80th percentiles of the same histogram, so the shape of the distribution is readable and not just its centre. These are descriptive rather than diagnostic, because size selection moves all three percentiles together, and the spread between them does not separate library types either. The comparison that does work is a csRNA library against its own matched input, rather than one row of the table against another. |
| `ModeReadLength` | Commonest read length in the same histogram, the shorter length winning a tie. This one is weak on its own, because a healthy csRNA library is nearly flat in length and the winning bin flips between replicates of the same arm. A mode sitting far from the median is best read as a prompt to open the histogram, rather than as a measurement in itself. |
| `ModeReadFraction` | Percentage of reads that sit at `ModeReadLength`, so the height of the tallest bin. Of the length columns this is the one that separates library types cleanly. Nothing gates on it, and it is not independent evidence, because it largely restates `PctNuclear` and `csFRiP`. |
| `TSSDistinctFrag` | Distinct fragments inside the final TSS set, where a fragment is a `(position, strand, length)` triple from the tag directory, so both ends define it. Depends on sequencing depth, so it is the denominator for the two below rather than something to compare across libraries. The rule excludes organelle chromosomes, as `csFRiP` does. |
| `TSSChao1` | Chao1 extrapolation of the fragments the library could yield at infinite depth, from its singleton and doubleton counts. It needs no subsampling, which is why the table reports it in place of a depth-matched count. |
| `TSSSaturation` | Good's coverage, `1 - singletons/reads`: the estimated chance the next read comes from a fragment already seen. This is the number to read when you are deciding whether more sequencing would pay. |
| `snRNA5pPct` | Share of the signal within 500 bp of an annotated snRNA 5′ end that sits within 5 bp of it, which requires a stranded `qc / snrnas`. Sm-class snRNAs are capped Pol II transcripts, so a working library puts nearly all of it on the 5′ end. Unlike `PhosEfficiency` it needs no matched input library. |
| `tRNA5pPct` | The same measurement at mature tRNA 5′ ends, which requires a stranded `qc / trnas`. RNase P leaves a 5′-monophosphate there, which becomes ligation-competent when the phosphatase step fails, so this moves opposite to `snRNA5pPct`. The table carries it as a cross-check on which chemistry failed rather than as a gate, since it does not separate good libraries from bad ones on its own. |
| `csFRiP` | Fraction of nuclear reads falling in csRNA-called TSSs. It should be above 0.9 in a good csRNA library. |
| `sFRiP` | Fraction of nuclear reads in input-called (small RNA) TSSs. |
| `PctNuclear` | Percentage of reads on nuclear chromosomes, so excluding `qc / organelle_chroms`. |
| `csEnrichment` | Ratio of csRNA FRiP to input FRiP, so how enriched the capped initiation signal is over background. |
| `sDepletion` | Inverse ratio of small RNA signal between csRNA and input libraries. |
| `PretRNAPct` | Pre-tRNA reads over pre-tRNA reads plus reads in TSSs, so contaminant per unit of signal rather than per library read. Requires `qc / trnas`. |
| `PhosEfficiency` | Ratio of pre-tRNA depletion in csRNA against input, measuring 5′-phosphate removal. |
| `miRNADepletion` | Ratio of miRNA depletion in csRNA against input, an independent phosphorylation-efficiency metric. Requires `qc / mirnas`. |
| `miRNAPct` | miRNA reads over miRNA reads plus reads in TSSs, the same basis as `PretRNAPct`, so the two contamination percentages are directly comparable. |
| `StrandBalance` | Fraction of mapped reads on the + strand, around 0.5 in csRNA libraries. Large deviations flag adapter contamination, strand bias in library prep, or pile-ups at a few highly expressed loci. Looser bounds in input libraries, where small-RNA biology is genuinely strand-skewed. |
| `TSSDetected` | Fraction of merged-set TSSs with at least one tag in this library. Low values flag undersequenced libraries. |
| `FinalTSSDetected` | Final tables only. Fraction of `tss.final.bed` TSSs with at least one tag in this library. |
| `NConsensusTSS`, `NFinalTSS`, `NFilteredTSS` | Final tables only. TSSs in `tss.consensus.bed`, TSSs retained in `tss.final.bed`, and TSSs the filters dropped. |
| `MinReplCorrSpearman` | Lowest Spearman correlation between this sample's TSS counts and any other csRNA replicate of the same `sample_name`. `NA` with only one replicate. Sharp drops, below 0.9 say, flag sample swaps and replicate dropouts. csRNA only. |
| `MinReplCorrPearson` | The same on log1p-transformed counts. |
| `Status` | `Ok` or `FAIL`, from the gates named in `StatusReason`. A gate the pipeline cannot evaluate, an `NA` csFRiP say, counts as a failure rather than a pass. csRNA only. |
| `StatusReason` | Comma-separated gates that tripped, empty when `Status` is `Ok`. |

## Reading these numbers

**The same column name covers different regions in the two files.** `qc_initial_*.txt`
measures `csRiP`, `csFRiP` and `csEnrichment` against the initial merged TSS set, and
`qc_final_*.txt` measures the same three against `tss.final.bed`. Both are correct, and
they are not the same measurement, so never join or compare the two files column by
column.

**The contamination percentages are signal-relative, not library-relative.**
`PretRNAPct` and `miRNAPct` are contaminant reads over contaminant plus reads-in-TSS, so
they answer "per read I care about, how much junk", not "what fraction of the library is
junk". The library-fraction version is smaller by roughly the csFRiP factor, which is
itself lower in worse libraries, so a contaminated library looks somewhat worse here
than it would on a library-fraction metric.

**A high duplicate rate is not automatically bad, and a low one is not automatically
good.** Duplication at matched depth is complexity convolved with concentration, since
reads pile up where the signal is, so a library that puts more of its reads into TSSs
duplicates more at the same depth. Read `TSSDistinctFrag`, `TSSChao1` and
`TSSSaturation` instead. The pipeline restricts all three to the called TSS set for
exactly this reason.

**Every per-class figure counts uniquely-mapping reads only.** All of them come from the
filtered BAM, which is past the MAPQ floor, so `PretRNAPct`, `miRNAPct` and any rRNA
figure measure the unique-mapping tail of those species rather than their abundance.
Assemblies commonly represent the rDNA arrays as many near-identical copies, which no
short read can place uniquely, so a suspiciously low rRNA figure is better read as
evidence that the rRNA sits in the discarded multi-mapping fraction. `BelowMapqReads`
says how large that fraction is, and `program / keep_below_mapq_sample` will hand you
the reads to identify.

**`csEnrichment` does not discriminate at small group sizes.** On two to three libraries
per condition it failed to separate a treatment that demonstrably changed the libraries.
`PretRNAPct`, `PhosEfficiency`, alignment survival and the useful-read fraction did
separate consistently, so power comparisons belong on those.

**Only the initial `Status` acts on anything.** `collect_consensus_tss` reads
`qc_initial_cs.txt`, so that table's gates, `qc / min_cs_frip` and `qc /
min_pct_nuclear`, are the ones that can change a run, and `filtering /
exclude_failed_from_consensus` is what makes them do so. The final table's `Status` is
advice with its own gates, because it measures a different thing, and nothing reads it,
so a FAIL there has already had no effect on the TSS set that table describes.

## Supporting QC files

| File | Description |
|------|-------------|
| `qc/replicate_correlation.txt` | Pairwise Spearman and Pearson correlations between csRNA replicates of one `sample_name`, for the initial TSS set (`Stage=Initial`) and the filtered set (`Stage=Final`). Empty when no `sample_name` has two or more replicates. |
| `qc/stats_initial_cs.txt`, `qc/stats_initial_in.txt` | Raw alignment and tag counts per library: `RawReads`, `TrimmedReads`, `AlignedReads`, `BelowMapqReads`, `FilteredOutReads`, `TotalReads`, the read-length summary, `OrganelleReads`, `PosReads`, `NegReads`. |
| `tss/{sample}.stats.txt` | HOMER's own report from TSS calling: cluster counts, the enrichment threshold it chose (`log2 fold vs. input`), and the true and false positive set sizes it chose from. The pipeline post-processes this file. It rewrites the promoter-distal and stable-transcript fractions to `na` when HOMER had no annotation and no RNA-seq to compute them from. HOMER otherwise reports its placeholders (100.00% and 0.00%) as measurements. Expect `na` in those two positions. |
| `qc/{sample}.trimming.txt` | bfqutils trimming summary. |
| `qc/{sample}.aln.txt` | samtools flagstat on the **filtered** BAM, so it reports 100% mapped for every library and says nothing about what the filter removed. `aln.raw.txt` is the file that does. |
| `qc/{sample}.aln.raw.txt` | Accounting of the unfiltered alignment. It counts the reads seen, then the unmapped, secondary and supplementary ones, then those with a primary alignment. For that last group it says how many fell below the MAPQ floor, how many failed the remaining filters, and how many passed. Ends with a MAPQ histogram and the filter settings in force. Counted in the pipe, since the unfiltered alignment never reaches disk. |
| `qc/{sample}.unmapped.sample.fastq.gz`, `qc/{sample}.belowmapq.sample.fastq.gz` | Samples of discarded reads. The pipeline writes these only when the relevant `program / keep_*_sample` key is above zero. See [configuration.md](configuration.md). |
| `qc/{sample}.aligner.log` | Aligner stderr. For STAR, the rule appends `Log.final.out` to it. |
