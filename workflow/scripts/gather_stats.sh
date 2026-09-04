#!/usr/bin/env bash
set -euo pipefail
printf "Sample\tRawReads\tTrimmedReads\tAlignedReads\tBelowMapqReads\tFilteredOutReads\tTotalReads\tP20ReadLength\tMedianReadLength\tP80ReadLength\tModeReadLength\tOrganelleReads\tFreq1A\tPosReads\tNegReads\n" > "$OUT"
for i in "$@"; do
    if [[ ! -f "$OUT_PREFIX/$i/tagInfo.txt" ]]; then
        echo "ERROR: tagInfo.txt not found for sample $i (${OUT_PREFIX}/${i}/tagInfo.txt)" >&2
        exit 1
    fi
    bg_pos="$BG_PREFIX/${i}.raw.pos.bedGraph.gz"
    bg_neg="$BG_PREFIX/${i}.raw.neg.bedGraph.gz"
    if [[ ! -f "$bg_pos" || ! -f "$bg_neg" ]]; then
        echo "ERROR: raw bedGraphs not found for sample $i ($bg_pos, $bg_neg)" >&2
        exit 1
    fi
    printf "%s\t" "$i" >> "$OUT"
    # Sum all "Processed N reads" lines (multiple lines when input had several files)
    awk '/^Processed/{gsub(",","",$0); s+=$2} END{print s+0}' "$LOG_PREFIX/${i}.trimming.txt" | tr -d '\n' >> "$OUT"
    printf "\t" >> "$OUT"
    awk '/^Output/{gsub(",","",$0); s+=$2} END{print s+0}' "$LOG_PREFIX/${i}.trimming.txt" | tr -d '\n' >> "$OUT"
    printf "\t" >> "$OUT"
    # Where the reads between TrimmedReads and TotalReads went, from aln_summary.awk,
    # which counted the unfiltered alignment stream: reads with a primary alignment, then
    # how many of those the MAPQ filter and the remaining filters (flags, alignment
    # length, mismatches) discard. The flagstat in ${i}.aln.txt cannot answer this because
    # it runs on the already-filtered BAM.
    aln_raw="$LOG_PREFIX/${i}.aln.raw.txt"
    if [[ ! -f "$aln_raw" ]]; then
        echo "ERROR: alignment summary not found for sample $i ($aln_raw)" >&2
        exit 1
    fi
    for key in PrimaryMapped BelowMapq FailedOtherFilters; do
        awk -F'\t' -v k="$key" '$1 == k {print $2; found = 1}
            END {if (!found) print "NA"}' "$aln_raw" | tr -d '\n' >> "$OUT"
        printf "\t" >> "$OUT"
    done
    grep "^genome=" "$OUT_PREFIX/$i/tagInfo.txt" | cut -f3 | tr -d '\n' >> "$OUT"
    printf "\t" >> "$OUT"
    # Read-length summary from HOMER's own length histogram, which is built from the
    # filtered BAM, so these describe the reads that survived into the tagdir rather
    # than everything the trimmer emitted, and the average HOMER prints in that file's
    # header is pulled toward the tail by a handful of long reads.
    #
    # The median is the stable number. Read the mode with care: a healthy csRNA library
    # has a nearly flat distribution, top bin around 5 to 6 per cent, so which length
    # wins flips between replicates of the same arm and means little on its own. It is
    # worth having because a library with a real discrete population is not flat at all
    # (22 per cent at 24 nt in an input library, 15 per cent at 30 nt in a phosphatase-
    # free arm), and in those the mode lands where the population is.
    #
    # P20 and P80 are here so the shape can be read rather than just the centre. The
    # useful quantity is the spread, not either percentile on its own: size selection
    # moves all three together, so a 40-70 nt input and a 20-70 nt csRNA library can
    # land on the same P80 while having quite different distributions. P80 minus the
    # median is what says the library is top-heavy, and a csRNA library run beside its
    # own input reliably carries the longer upper tail.
    #
    # HOMER writes the histogram in ascending length order, which the single pass below
    # relies on for all three percentiles.
    lendist="$OUT_PREFIX/$i/tagLengthDistribution.txt"
    if [[ -f "$lendist" ]]; then
        awk -F'\t' '
            $1 ~ /^[0-9]+$/ { n++; len[n] = $1 + 0; frac[n] = $2 + 0; total += $2 + 0 }
            END {
                if (n == 0 || total <= 0) { printf "NA\tNA\tNA\tNA"; exit }
                p20 = "NA"; med = "NA"; p80 = "NA"; mode = "NA"; best = -1; cum = 0
                for (k = 1; k <= n; k++) {
                    cum += frac[k]
                    if (p20 == "NA" && cum >= 0.2 * total) p20 = len[k]
                    if (med == "NA" && cum >= 0.5 * total) med = len[k]
                    if (p80 == "NA" && cum >= 0.8 * total) p80 = len[k]
                    # > rather than >=, so a tie keeps the shorter length
                    if (frac[k] > best) { best = frac[k]; mode = len[k] }
                }
                printf "%s\t%s\t%s\t%s", p20, med, p80, mode
            }' "$lendist" >> "$OUT"
    else
        printf "NA\tNA\tNA\tNA" >> "$OUT"
    fi
    printf "\t" >> "$OUT"
    # Summed and rounded inside awk, not in the shell: HOMER writes tag counts as floats,
    # and a fractional total (multi-mapper weighting, bam2td -keepAll) comes back from awk
    # in %.6g scientific notation, which bash integer arithmetic rejects outright with
    # "syntax error: invalid arithmetic operator".
    organelle_total=$(awk -F'\t' -v chroms="$ORGANELLE_CHROMS" '
        BEGIN { n = split(chroms, a, /[[:space:]]+/)
                for (k = 1; k <= n; k++) if (a[k] != "") want[a[k]] = 1 }
        ($1 in want) { sum += $3 }
        END { printf "%.0f", sum + 0 }' "$OUT_PREFIX/$i/tagInfo.txt")
    printf "%s\t" "$organelle_total" >> "$OUT"
    awk '$1 == 0 {print $2}' "$OUT_PREFIX/$i/tagFreq.txt" | tr -d '\n' >> "$OUT"
    printf "\t" >> "$OUT"
    # Sum raw read counts on each strand from the bedGraphs. makeUCSCfile -raw writes
    # positive values on both strands, so no sign handling is needed here; the negative
    # values are introduced later, in generate_normalized_bw.R, for display only. Each
    # bedGraph interval contributes value * (end - start) reads.
    gunzip -c "$bg_pos" | awk 'NF >= 4 && $2 ~ /^[0-9]+$/ {sum += $4 * ($3 - $2)} END {printf "%.0f", sum+0}' >> "$OUT"
    printf "\t" >> "$OUT"
    gunzip -c "$bg_neg" | awk 'NF >= 4 && $2 ~ /^[0-9]+$/ {sum += $4 * ($3 - $2)} END {printf "%.0f", sum+0}' >> "$OUT"
    printf "\n" >> "$OUT"
done
