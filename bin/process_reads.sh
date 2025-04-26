#!/usr/bin/env bash
#
# Author: Benjamin Jean-Marie Tremblay (benjamin.tremblay@tsl.ac.uk)
# Date created: 23 April 2025
# Date modified: 24 April 2025
#

set -eo pipefail

WD=`pwd`

INPUT=
INDEX=
ADP=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
MAX_LEN=70
BWA_THREADS=4
SAM_THREADS=4
BWA_MAPQ=30
OUT_QC="./qc"
OUT_TRIM="./trim"
OUT_BAM="./bam"
OUT_TAGDIR="./tagdir"
SAMPLE_SHEET="./samples-processed.txt"
FORCE=

alias echo='echo "["$(date +"%Y-%m-%d %T")"]" '

help() {
    echo "csRNA-seq pre-processing: Process raw reads & gather QC data"
    echo "Benjamin Jean-Marie Tremblay (2025)"
    echo
    echo "Requirements: homer (+tair10), samtools, bfqutils, bwa"
    echo
    echo "Options:"
    echo "-i   TSV with 3 columns: 1) [sample name]_[cs/in][rep #],"
    echo "     2) path to read 1 FASTQ, 3) path to read 2 FASTQ (if PE reads)."
    echo "-o   Output sample sheet for the next step. (default=${SAMPLE_SHEET})"
    echo "-m   Max read length. (default=${MAX_LEN})"
    echo "-a   Adapter sequence for trimming SE reads. (default=${ADP})"
    echo "-I   Path to genome in bwa index folder."
    echo "-j   Number of threads for bwa-aln. (default=${BWA_THREADS})"
    echo "-J   Number of threads for samtools-sort. (default=${SAM_THREADS})"
    echo "-Q   bwa-aln mapping score cutoff for aligned reads. (default=${BWA_MAPQ})"
    echo "-q   Output dir for QC files. (default=${OUT_QC})"
    echo "-T   Output dir for trimmed FASTQ files. (default=${OUT_TRIM})"
    echo "-b   Output dir for BAM files. (default=${OUT_BAM})"
    echo "-t   Output dir for tag directories. (default=${OUT_TAGDIR})"
    echo "-f   Force re-run of alignment and/or tag directory creation even"
    echo "     if existing files/folders are present."
    echo "-h   Show this message."
    echo
    echo "If using containers, set the following variables in your env:"
    echo '$HOMER, $SAMTOOLS, $BFQUTILS, $BWA'
    echo '(e.g. export HOMER="singularity exec /path/to/containers/homer.img")'
    exit 1
}

while getopts "i:o:a:m:I:j:Q:q:T:b:t:fh" opt ; do
    case $opt in 
        i) INPUT=$OPTARG;;
        o) SAMPLE_SHEET=$OPTARG;;
        a) ADP=$OPTARG;;
        m) MAX_LEN=$OPTARG;;
        I) INDEX=$OPTARG;;
        j) BWA_THREADS=$OPTARG;;
        Q) BWA_MAPQ=$OPTARG;;
        q) OUT_QC=$OPTARG;;
        T) OUT_TRIM=$OPTARG;;
        b) OUT_BAM=$OPTARG;;
        t) OUT_TAGDIR=$OPTARG;;
        f) FORCE="force";;
        h) help;;
        ?) help;;
    esac
done

[ -z $INPUT ] && { echo "Error: Missing -i"; exit 1; }
[ -z $INDEX ] && { echo "Error: Missing -I"; exit 1; }

cat $INPUT | tr -d '\r' > ./tmp_input

cleanup() {
    cd $WD
    rm -f ./tmp_input
}
trap cleanup EXIT

touch ./tmp_err

printErr() {
    cd $WD
    echo "See error log: `realpath ./tmp_err`"
    echo "Printing last lines ..."
    echo "--------------"
    tail -n 5 ./tmp_err
    echo "--------------"
}
trap printErr ERR

mkdir -p ${OUT_QC}
mkdir -p ${OUT_TRIM}
mkdir -p ${OUT_BAM}
mkdir -p ${OUT_TAGDIR}

rm -f ${SAMPLE_SHEET}
touch ${SAMPLE_SHEET}

echo "Processing `grep . ./tmp_input | wc -l` samples ..."

while IFS=$'\t' read -r SAMPLE R1 R2 || [ $SAMPLE ] ; do

    [ -z $SAMPLE ] && continue

    echo "    Working on sample: $SAMPLE"
    [ -z $R1 ] && { echo "Error: missing R1 fastq for sample $SAMPLE"; exit 1; }
    [ -z $R2 ] && echo "        Found no R2 fastq, assuming SE."

    TOTAL_READS=0
    FINAL_READS=0

    # TODO: Parallelize trimming with xargs across multiple samples simultaneously!
    if [ -z $R2 ] ; then
        if [ ! -f ${OUT_TRIM}/${SAMPLE}.trim.fq.gz ] || [ ! -z $FORCE ] ; then
            echo "        Trimming SE reads ..."
            gunzip -c $R1 | $BFQUTILS bfqtrimse -m $MAX_LEN -a $ADP - \
                2> ${OUT_QC}/${SAMPLE}.trimming.txt \
            | $BFQUTILS bfqstats -O -N -z \
                -l ${OUT_QC}/${SAMPLE}.trim.lengths.txt \
                -g ${OUT_QC}/${SAMPLE}.trim.gc.txt \
                -q ${OUT_QC}/${SAMPLE}.trim.qual.txt \
                -Q ${OUT_QC}/${SAMPLE}.trim.ppqual.txt \
                -b ${OUT_QC}/${SAMPLE}.trim.ppbase.txt \
                -k ${OUT_QC}/${SAMPLE}.trim.kmer.txt \
                -o ${OUT_QC}/${SAMPLE}.trim.summary.txt \
                > ${OUT_TRIM}/${SAMPLE}.trim.fq.gz
        else
            echo "        Skipping trimming SE reads"
        fi
    else
        if [ ! -f ${OUT_TRIM}/${SAMPLE}.trim.fq.gz ] || [ ! -z $FORCE ] ; then
            echo "        Merging PE reads ..."
            $BFQUTILS bfqmerge -m $MAX_LEN <(gunzip -c $R1) <(gunzip -c $R2) \
                2> ${OUT_QC}/${SAMPLE}.trimming.txt \
            | $BFQUTILS bfqstats -O -N -z \
                -l ${OUT_QC}/${SAMPLE}.trim.lengths.txt \
                -g ${OUT_QC}/${SAMPLE}.trim.gc.txt \
                -q ${OUT_QC}/${SAMPLE}.trim.qual.txt \
                -Q ${OUT_QC}/${SAMPLE}.trim.ppqual.txt \
                -b ${OUT_QC}/${SAMPLE}.trim.ppbase.txt \
                -k ${OUT_QC}/${SAMPLE}.trim.kmer.txt \
                -o ${OUT_QC}/${SAMPLE}.trim.summary.txt \
                > ${OUT_TRIM}/${SAMPLE}.trim.fq.gz
        else
            echo "        Skipping merging PE reads"
        fi
    fi

    TOTAL_READS=`head -n 1 ${OUT_QC}/${SAMPLE}.trimming.txt | awk '{print $2}'`
    FINAL_READS=`tail -n 1 ${OUT_QC}/${SAMPLE}.trimming.txt | awk '{print $2}'`

    if [ ! -f ${OUT_BAM}/${SAMPLE}.sai ] || [ ! -z $FORCE ] ; then
        echo "        Aligning reads ..."
        $BWA bwa aln -o 0 -t ${BWA_THREADS} ${INDEX} \
            ${OUT_TRIM}/${SAMPLE}.trim.fq.gz 2> tmp_err > ${OUT_BAM}/${SAMPLE}.sai
    else
        echo "        Skipping aligning reads"
    fi
    if [ ! -f ${OUT_BAM}/${SAMPLE}.sort.bam ] || [ ! -z $FORCE ] ; then
        echo "        Filtering and sorting alignments ..."
        $BWA bwa samse ${INDEX} ${OUT_BAM}/${SAMPLE}.sai ${INDEX} \
            ${OUT_TRIM}/${SAMPLE}.trim.fq.gz 2> tmp_err \
        | $SAMTOOLS samtools view -F 2308 -q ${BWA_MAPQ} -b - \
        | $SAMTOOLS samtools sort --verbosity 0 -@ ${SAM_THREADS} \
            -o ${OUT_BAM}/${SAMPLE}.sort.bam -
        echo "        Getting alignment stats ..."
        $SAMTOOLS samtools index ${OUT_BAM}/${SAMPLE}.sort.bam 2> tmp_err
        $SAMTOOLS samtools flagstat ${OUT_BAM}/${SAMPLE}.sort.bam 2> tmp_err \
            > ${OUT_QC}/${SAMPLE}.aln.txt
    else
        echo "        Skipping filtering and sorting alignments"
    fi

    if [ ! -d ${OUT_TAGDIR}/${SAMPLE} ] || [ ! -z $FORCE ] ; then
        echo "        Creating tag directory ..."
        $SAMTOOLS samtools view -h ${OUT_BAM}/${SAMPLE}.sort.bam 2> tmp_err \
            > ${OUT_BAM}/${SAMPLE}.sort.sam
        $HOMER makeTagDirectory ${OUT_TAGDIR}/${SAMPLE} \
            ${OUT_BAM}/${SAMPLE}.sort.sam -genome tair10 -checkGC 2> tmp_err
        rm ${OUT_BAM}/${SAMPLE}.sort.sam
    else
        echo "        Skipping creating tag directory"
    fi

    if [ `echo ${SAMPLE} | grep "_cs[0-9][0-9]*"` ] ; then
        SAMPLE_NAME=`echo ${SAMPLE} | sed 's/_cs[0-9][0-9]*//g'`
        SAMPLE_REP=`echo ${SAMPLE} | grep -o "_cs[0-9][0-9]*" | grep -o "[0-9][0-9]*"`
        printf "${SAMPLE_NAME}\t${SAMPLE_REP}\t${OUT_TAGDIR}\n" >> ${SAMPLE_SHEET}
    fi

    MAPPED_READS=`grep "primary mapped" ${OUT_QC}/$SAMPLE}.aln.txt | awk '{print $1}'`

    echo "    >>> Finished sample ${SAMPLE}"
    if [ -z $R2 ] ; then
        echo "        Starting read count: ${TOTAL_READS}"
    else 
        echo "        Starting read pairs: ${TOTAL_READS}"
    fi
    echo "        Final read count: ${FINAL_READS}"
    echo "        Aligned reads: ${MAPPED_READS}"

done < ./tmp_input

echo "All done."

rm ./tmp_err

