#!/usr/bin/env bash
set -euo pipefail
printf "Sample\tTotalReads\tOrganelleReads\tFreq1A\n" > "$OUT"
for i in "$@"; do
    if [[ ! -f "$OUT_PREFIX/$i/tagInfo.txt" ]]; then
        echo "ERROR: tagInfo.txt not found for sample $i (${OUT_PREFIX}/${i}/tagInfo.txt)" >&2
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
    printf "\n" >> "$OUT"
done
