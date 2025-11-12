printf "Sample\tTotalReads\tPtReads\tMtReads\tFreq1A\n" > $OUT
for i in "$@"; do
    printf "$i\t" >> $OUT
    grep "^genome=" $OUT_PREFIX/$i/tagInfo.txt | cut -f3 | tr -d '\n' >> $OUT
    printf "\t" >> $OUT
    if [ `grep "^Pt" $OUT_PREFIX/$i/tagInfo.txt | wc -l` -gt 0 ] ; then
        grep "^Pt" $OUT_PREFIX/$i/tagInfo.txt | cut -f3 | tr -d '\n' >> $OUT
    else
        printf "0" >> $OUT
    fi
    printf "\t" >> $OUT
    if [ `grep "^Mt" $OUT_PREFIX/$i/tagInfo.txt | wc -l` -gt 0 ] ; then
        grep "^Mt" $OUT_PREFIX/$i/tagInfo.txt | cut -f3 | tr -d '\n' >> $OUT
    else
        printf "0" >> $OUT
    fi
    printf "\t" >> $OUT
    grep "^0" $OUT_PREFIX/$i/tagFreq.txt | cut -f2 | tr -d '\n' >> $OUT
    printf "\n" >> $OUT
done
