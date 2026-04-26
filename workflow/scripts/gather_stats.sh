printf "Sample\tTotalReads\tOrganelleReads\tFreq1A\n" > $OUT
for i in "$@"; do
    printf "$i\t" >> $OUT
    grep "^genome=" $OUT_PREFIX/$i/tagInfo.txt | cut -f3 | tr -d '\n' >> $OUT
    printf "\t" >> $OUT
    organelle_total=0
    for chrom in $ORGANELLE_CHROMS; do
        n=$(awk -v c="$chrom" '$0 ~ "^"c"\t" {sum += $3} END {print sum+0}' $OUT_PREFIX/$i/tagInfo.txt)
        organelle_total=$((organelle_total + n))
    done
    printf "%d" $organelle_total >> $OUT
    printf "\t" >> $OUT
    grep "^0" $OUT_PREFIX/$i/tagFreq.txt | cut -f2 | tr -d '\n' >> $OUT
    printf "\n" >> $OUT
done
