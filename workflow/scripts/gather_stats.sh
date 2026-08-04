#!/usr/bin/env bash
set -euo pipefail
printf "Sample\tRawReads\tTrimmedReads\tTotalReads\tOrganelleReads\tFreq1A\tPosReads\tNegReads\n" > "$OUT"
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
    grep "^genome=" "$OUT_PREFIX/$i/tagInfo.txt" | cut -f3 | tr -d '\n' >> "$OUT"
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
