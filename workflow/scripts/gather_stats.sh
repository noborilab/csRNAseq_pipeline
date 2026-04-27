#!/usr/bin/env bash
set -euo pipefail
printf "Sample\tTotalReads\tOrganelleReads\tFreq1A\tPosReads\tNegReads\n" > "$OUT"
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
    grep "^genome=" "$OUT_PREFIX/$i/tagInfo.txt" | cut -f3 | tr -d '\n' >> "$OUT"
    printf "\t" >> "$OUT"
    organelle_total=0
    for chrom in $ORGANELLE_CHROMS; do
        n=$(awk -v c="$chrom" 'index($0, c "\t") == 1 {sum += $3} END {print sum+0}' "$OUT_PREFIX/$i/tagInfo.txt")
        organelle_total=$((organelle_total + n))
    done
    printf "%d\t" $organelle_total >> "$OUT"
    awk '$1 == 0 {print $2}' "$OUT_PREFIX/$i/tagFreq.txt" | tr -d '\n' >> "$OUT"
    printf "\t" >> "$OUT"
    # Sum raw read counts on each strand from the bedGraphs (negative strand
    # values are stored as negative numbers by HOMER, so negate to recover the
    # positive read count). Each bedGraph interval is value * (end - start) reads.
    gunzip -c "$bg_pos" | awk 'NF >= 4 && $2 ~ /^[0-9]+$/ {sum += $4 * ($3 - $2)} END {printf "%.0f", sum+0}' >> "$OUT"
    printf "\t" >> "$OUT"
    gunzip -c "$bg_neg" | awk 'NF >= 4 && $2 ~ /^[0-9]+$/ {sum += $4 * ($3 - $2)} END {printf "%.0f", sum+0}' >> "$OUT"
    printf "\n" >> "$OUT"
done
