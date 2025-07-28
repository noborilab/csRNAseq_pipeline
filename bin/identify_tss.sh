#!/usr/bin/env bash
#
# Author: Benjamin Jean-Marie Tremblay (benjamin.tremblay@tsl.ac.uk)
# Date created: 10 March 2025
# Date modified: 28 July 2025
#

set -eo pipefail

WD=`pwd`

INPUT=
OUT_BG="./bedGraphs"
OUT_TSS="./tss"
OUT_BED="./tss_bed"
SAMPLE_SHEET="./samples-tss.txt"
OUT_FINAL="."
SKIP_TSS=
SKIP_BG=
FORCE=
USE_REL=

help() {
    echo "csRNA-seq post-processing part 1: Identify TSSs, quantify them, and make bedGraphs"
    echo "Benjamin Jean-Marie Tremblay (2025)"
    echo
    echo "Requirements: homer (+tair10), bedtools"
    echo
    echo "Options:"
    echo "-i   TSV with 3 columns: 1) sample name, 2) rep #, 3) tagdir path."
    echo "-o   Output sample sheet for finalizing TSSs script. (default=${SAMPLE_SHEET})"
    echo "-g   Output dir for bedGraphs (default=${OUT_BG})"
    echo "-t   Output dir for HOMER TSS results (default=${OUT_TSS})"
    echo "-b   Output dir for TSS BED files (default=${OUT_BED})"
    echo "-O   Output dir for final TSSs, quantification data and stats (default=${OUT_FINAL})"
    echo "-s   Skip TSS identification"
    echo "-b   Skip bedGraph creation"
    echo "-f   Force re-run of TSS ID & bedGraphs steps for which output files already exist."
    echo "-r   Use relative paths for output sample sheet."
    echo "-h   Show this message"
    echo
    echo "If using containers, set the following variables in your env:"
    echo '$HOMER, $BEDTOOLS'
    echo '(e.g. export HOMER="singularity exec /path/to/containers/homer.img")'
    exit 1
}

while getopts "i:o:g:t:b:O:sbfrh" opt; do
    case $opt in
        i) INPUT=$OPTARG;;
        o) SAMPLE_SHEET=$OPTARG;;
        g) OUT_BG=$OPTARG;;
        t) OUT_TSS=$OPTARG;;
        b) OUT_BED=$OPTARG;;
        O) OUT_FINAL=$OPTARG;;
        s) SKIP_TSS="skip";;
        b) SKIP_BG="skip";;
        f) FORCE="force";;
        r) USE_REL="userel";;
        h) help;;
        ?) help;;
    esac
done

[ -z $INPUT ] && { echo "Error: Missing -i"; exit 1; }

cat $INPUT | tr -d '\r' > ./tmp_input

cleanup() {
    cd $WD
    rm -f ./tmp_input
    rm -f ./tmp_tagdirs_cs
    rm -f ./tmp_tagdirs_in
}
trap cleanup EXIT

touch ./tmp_err

printErr() {
    cd $WD
    echo "ERROR: Stopping"
    echo "See error log: `realpath ./tmp_err`"
    echo "Printing last lines ..."
    echo "--------------"
    tail -n 5 ./tmp_err
    echo "--------------"
}
trap printErr ERR

mkdir -p $OUT_BG
mkdir -p $OUT_TSS
mkdir -p $OUT_BED

rm -f ${OUT_BED}/all_cs.tss.bed
touch ${OUT_BED}/all_cs.tss.bed

rm -f ${OUT_BED}/all_in.tss.bed
touch ${OUT_BED}/all_in.tss.bed

rm -f ./tmp_tagdirs_cs
touch ./tmp_tagdirs_cs

rm -f ./tmp_tagdirs_in
touch ./tmp_tagdirs_in

rm -f ${SAMPLE_SHEET}
touch ${SAMPLE_SHEET}

printf "Sample\tTotalReads\tPtReads\tMtReads\tFreq1A\n" > ${OUT_FINAL}/stats_cs.txt
printf "Sample\tTotalReads\tPtReads\tMtReads\tFreq1A\n" > ${OUT_FINAL}/stats_in.txt

echo "Identifying TSSs & making bedGraphs ..."

while IFS=$'\t' read -r SAMPLE REP TAGDIRS || [ $SAMPLE ] ; do
    [ -z $SAMPLE ] && continue
    [ -z $REP ] && { echo "Error: missing rep # for sample $SAMPLE"; exit 1; }
    [ -z $TAGDIRS ] && { echo "Error: missing tagdirs path for sample $SAMPLE"; exit 1; }
    echo "    Working on sample: $SAMPLE r${REP}"

    if [ -z $SKIP_TSS ] ; then
        if [ ! -f ${OUT_BED}/${SAMPLE}_cs${REP}.tss.bed ] || [ ! -z $FORCE ] ; then
            $HOMER findcsRNATSS.pl ${TAGDIRS}/${SAMPLE}_cs${REP} \
                -o ${OUT_TSS}/${SAMPLE}_cs${REP} \
                -i ${TAGDIRS}/${SAMPLE}_in${REP} \
                -genome tair10 \
                -size 100 \
                -ntagThreshold 5 \
                -includeSingleExons 2> tmp_err
            $HOMER pos2bed.pl -o ${OUT_BED}/${SAMPLE}_cs${REP}.tss.bed \
                ${OUT_TSS}/${SAMPLE}_cs${REP}.tss.txt 2> tmp_err
        fi
    fi
    cat ${OUT_BED}/${SAMPLE}_cs${REP}.tss.bed >> ${OUT_BED}/all_cs.tss.bed

    if [ -z $SKIP_TSS ] ; then
        if [ ! -f ${OUT_BED}/${SAMPLE}_in${REP}.tss.bed ] || [ ! -z $FORCE ] ; then
            $HOMER findcsRNATSS.pl ${TAGDIRS}/${SAMPLE}_in${REP} \
                -o ${OUT_TSS}/${SAMPLE}_in${REP} \
                -i ${TAGDIRS}/${SAMPLE}_cs${REP} \
                -genome tair10 \
                -size 100 \
                -ntagThreshold 5 \
                -includeSingleExons 2> tmp_err
            $HOMER pos2bed.pl -o ${OUT_BED}/${SAMPLE}_in${REP}.tss.bed \
                ${OUT_TSS}/${SAMPLE}_in${REP}.tss.txt 2> tmp_err
        fi
    fi
    cat ${OUT_BED}/${SAMPLE}_in${REP}.tss.bed >> ${OUT_BED}/all_in.tss.bed

    NTSS_CS=`grep "Valid TSS clusters" ${OUT_TSS}/${SAMPLE}_cs${REP}.stats.txt | cut -f2`
    NTSS_IN=`grep "Valid TSS clusters" ${OUT_TSS}/${SAMPLE}_in${REP}.stats.txt | cut -f2`

    echo "        Found ${NTSS_CS} csRNA TSS and ${NTSS_IN} input TSS."


    TD_CS=${TAGDIRS}/${SAMPLE}_cs${REP}
    TD_IN=${TAGDIRS}/${SAMPLE}_in${REP}
    echo $TD_CS >> ./tmp_tagdirs_cs
    echo $TD_IN >> ./tmp_tagdirs_in

    printf "${SAMPLE}_cs${REP}\t" >> ${OUT_FINAL}/stats_cs.txt
    grep "^genome=" ${TD_CS}/tagInfo.txt | cut -f3 | tr -d '\n' >> ${OUT_FINAL}/stats_cs.txt
    printf "\t" >> ${OUT_FINAL}/stats_cs.txt
    grep "^Pt" ${TD_CS}/tagInfo.txt | cut -f3 | tr -d '\n' >> ${OUT_FINAL}/stats_cs.txt
    printf "\t" >> ${OUT_FINAL}/stats_cs.txt
    grep "^Mt" ${TD_CS}/tagInfo.txt | cut -f3 | tr -d '\n' >> ${OUT_FINAL}/stats_cs.txt
    printf "\t" >> ${OUT_FINAL}/stats_cs.txt
    grep "^0" ${TD_CS}/tagFreq.txt | cut -f2 | tr -d '\n' >> ${OUT_FINAL}/stats_cs.txt
    printf "\n" >> ${OUT_FINAL}/stats_cs.txt

    if [ -z $SKIP_BG ] ; then
        if [ ! -f ${OUT_BG}/${SAMPLE}_cs${REP}.raw.pos.bedGraph.gz ] || [ ! -z $FORCE ] ; then
            $HOMER makeUCSCfile ${TD_CS} -style tss -strand + -raw \
                -o ${OUT_BG}/${SAMPLE}_cs${REP}.raw.pos.bedGraph 2> tmp_err
        fi
        if [ ! -f ${OUT_BG}/${SAMPLE}_cs${REP}.raw.neg.bedGraph.gz ] || [ ! -z $FORCE ] ; then
            $HOMER makeUCSCfile ${TD_CS} -style tss -strand - -raw \
                -o ${OUT_BG}/${SAMPLE}_cs${REP}.raw.neg.bedGraph 2> tmp_err
        fi
    fi

    printf "${SAMPLE}_in${REP}\t" >> ${OUT_FINAL}/stats_in.txt
    grep "^genome=" ${TD_IN}/tagInfo.txt | cut -f3 | tr -d '\n' >> ${OUT_FINAL}/stats_in.txt
    printf "\t" >> ${OUT_FINAL}/stats_in.txt
    grep "^Pt" ${TD_IN}/tagInfo.txt | cut -f3 | tr -d '\n' >> ${OUT_FINAL}/stats_in.txt
    printf "\t" >> ${OUT_FINAL}/stats_in.txt
    grep "^Mt" ${TD_IN}/tagInfo.txt | cut -f3 | tr -d '\n' >> ${OUT_FINAL}/stats_in.txt
    printf "\t" >> ${OUT_FINAL}/stats_in.txt
    grep "^0" ${TD_IN}/tagFreq.txt | cut -f2 | tr -d '\n' >> ${OUT_FINAL}/stats_in.txt
    printf "\n" >> ${OUT_FINAL}/stats_in.txt

    if [ -z $SKIP_BG ] ; then
        if [ ! -f ${OUT_BG}/${SAMPLE}_in${REP}.raw.pos.bedGraph.gz ] || [ ! -z $FORCE ] ; then
            $HOMER makeUCSCfile ${TD_IN} -style tss -strand + -raw \
                -o ${OUT_BG}/${SAMPLE}_in${REP}.raw.pos.bedGraph 2> tmp_err
        fi
        if [ ! -f ${OUT_BG}/${SAMPLE}_in${REP}.raw.neg.bedGraph.gz ] || [ ! -z $FORCE ] ; then
            $HOMER makeUCSCfile ${TD_IN} -style tss -strand - -raw \
                -o ${OUT_BG}/${SAMPLE}_in${REP}.raw.neg.bedGraph 2> tmp_err
        fi
    fi

    if [ -z ${USE_REL} ] ; then 
        TAGDIRS_O=`realpath ${TAGDIRS}`
        OUT_BG_O=`realpath ${OUT_BG}`
        printf "${SAMPLE}\t${REP}\t${TAGDIRS_O}\t${OUT_BG_O}\n" >> ${SAMPLE_SHEET}
    else
        TAGDIRS_O=`realpath --relative-to=. ${TAGDIRS}`
        OUT_BG_O=`realpath --relative-to=. ${OUT_BG}`
        printf "${SAMPLE}\t${REP}\t./${TAGDIRS_O}\t./${OUT_BG_O}\n" >> ${SAMPLE_SHEET}
    fi

done < ./tmp_input

TAGDIRS_CS=(`cat ./tmp_tagdirs_cs`)
TAGDIRS_IN=(`cat ./tmp_tagdirs_in`)

echo "Merging TSS ..."
$BEDTOOLS bedtools sort -i ${OUT_BED}/all_cs.tss.bed \
    > ${OUT_FINAL}/all_cs.tss.sort.bed
$BEDTOOLS bedtools sort -i ${OUT_BED}/all_in.tss.bed \
    > ${OUT_FINAL}/all_in.tss.sort.bed

$BEDTOOLS bedtools merge -s -c 6 -o distinct -i ${OUT_FINAL}/all_cs.tss.sort.bed \
| $BEDTOOLS bedtools sort -i stdin \
| awk 'BEGIN{IFS="\t";OFS="\t"}{if ($1!="Mt" && $1!="Pt") {print $1,$2,$3,"TSS_" NR,"0",$4}}' \
    > ${OUT_FINAL}/all_cs.tss_merged.bed

$BEDTOOLS bedtools merge -s -c 6 -o distinct -i ${OUT_FINAL}/all_in.tss.sort.bed \
| $BEDTOOLS bedtools sort -i stdin \
| awk 'BEGIN{IFS="\t";OFS="\t"}{if ($1!="Mt" && $1!="Pt") {print $1,$2,$3,"TSS_" NR,"0",$4}}' \
    > ${OUT_FINAL}/all_in.tss_merged.bed

NTSS_CS_FINAL=`cat ${OUT_FINAL}/all_cs.tss_merged.bed | wc -l`
NTSS_IN_FINAL=`cat ${OUT_FINAL}/all_in.tss_merged.bed | wc -l`

echo "    Final csRNA TSS count: ${NTSS_CS_FINAL}"
echo "    Final input TSS count: ${NTSS_IN_FINAL}"

echo "Quantifying TSS ..."
$HOMER annotatePeaks.pl ${OUT_FINAL}/all_cs.tss_merged.bed \
    tair10 -strand + -fragLength 1 -raw -d ${TAGDIRS_CS[@]} \
    > ${OUT_FINAL}/all_cs.tss_merged_quant_cs.txt 2> ./tmp_err

$HOMER annotatePeaks.pl ${OUT_FINAL}/all_cs.tss_merged.bed \
    tair10 -strand + -fragLength 1 -raw -d ${TAGDIRS_IN[@]} \
    > ${OUT_FINAL}/all_cs.tss_merged_quant_in.txt 2> ./tmp_err

$HOMER annotatePeaks.pl ${OUT_FINAL}/all_in.tss_merged.bed \
    tair10 -strand + -fragLength 1 -raw -d ${TAGDIRS_CS[@]} \
    > ${OUT_FINAL}/all_in.tss_merged_quant_cs.txt 2> ./tmp_err

$HOMER annotatePeaks.pl ${OUT_FINAL}/all_in.tss_merged.bed \
    tair10 -strand + -fragLength 1 -raw -d ${TAGDIRS_IN[@]} \
    > ${OUT_FINAL}/all_in.tss_merged_quant_in.txt 2> ./tmp_err

echo "All done."

rm ./tmp_err

