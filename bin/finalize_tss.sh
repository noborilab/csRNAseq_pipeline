#!/usr/bin/env bash
#
# Author: Benjamin Jean-Marie Tremblay (benjamin.tremblay@tsl.ac.uk)
# Date created: 11 March 2025
# Date modified: 25 April 2025
#

set -e

WD=`pwd`

INPUT=
IN_TSS="./tss"
IN_BED="./bed"
OUT_BW="./bw"
OUT_FINAL="."
FORCE_QUANT=

help() {
    echo "csRNA-seq post-processing part 3: consolidate data from high quality samples"
    echo "Benjamin Jean-Marie Tremblay (2025)"
    echo
    echo "Requirements: homer (+tair10), bedtools, R (+rtracklayer, +readr, +edgeR)"
    echo
    echo "Options:"
    echo "-i   TSV with 4 columns: 1) sample name, 2) rep #, 3) tagdirs path, 4) bedGraphs path."
    echo "-t   Dir containing for per-sample TSSs (default=${IN_TSS})"
    echo "-b   Dir containing pre-sample TSS BED files (default=${IN_BED})"
    echo "-o   Output dir for final TSSs, and raw+normalized quantification tables (default=${OUT_FINAL})"
    echo "-w   Output dir for final bigWigs (default=${OUT_BW})"
    echo "-f   Force re-run of quantification even if the output file already exist."
    echo "-h   Show this message"
    echo
    echo "If using containers, set the following variables in your env:"
    echo '$HOMER, $BEDTOOLS, $RSCRIPT'
    echo '(e.g. export HOMER="singularity exec /path/to/containers/homer.img")'
    exit 1
}

while getopts "i:t:b:o:w:fh" opt; do
    case $opt in
        i) INPUT=$OPTARG;;
        t) IN_TSS=$OPTARG;;
        b) IN_BED=$OPTARG;;
        o) OUT_FINAL=$OPTARG;;
        w) OUT_BW=$OPTARG;;
        f) FORCE_QUANT="force";;
        h) help;;
        ?) help;;
    esac
done

[ -z $INPUT ] && { echo "Error: Missing -i"; exit 1; }

cat $INPUT | tr -d '\r' > ./tmp_input

rm -f ./tmp_tagdirs
touch ./tmp_tagdirs

cleanup() {
    cd $WD
    rm -f ./tmp_input
    rm -f ./tmp_tagdirs
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

mkdir -p $OUT_BW
mkdir -p $OUT_FINAL

while IFS=$'\t' read -r SAMPLE REP TAGDIRS BEDGRAPHS || [ $SAMPLE ] ; do
    [ -z $SAMPLE ] && continue
    [ -z $REP ] && { echo "Error: missing rep # for sample ${SAMPLE}_cs${REP}"; exit 1; }
    [ -z $TAGDIRS ] && { echo "Error: missing tagdir path for sample ${SAMPLE}_cs${REP}"; exit 1; }
    [ -z $BEDGRAPHS ] && { echo "Error: missing bedGraph path for sample ${SAMPLE}_cs${REP}"; exit 1; }

    echo ${TAGDIRS}/${SAMPLE}_cs${REP} >> ./tmp_tagdirs

done < ./tmp_input

echo "Finding consensus TSSs across reps and merging ..."

$RSCRIPT Rscript -e "
library(rtracklayer)
input <- read.table(trimws(\"$INPUT\"), stringsAsFactors=FALSE)
samples <- unique(input[[1]])
tss_all <- vector('list', length(samples))
for (i in seq_along(samples)) {
    input_i <- input[input[[1]] %in% samples[i], ]
    samples_i <- paste0(\"${IN_BED}\", '/', input_i[[1]], '_cs', input_i[[2]], '.tss.bed')
    tss_i <- import(samples_i[1])
    for (j in seq_along(samples_i[-1])) {
        tss_j <- import(samples_i[1 + j])
        tss_i <- tss_i[overlapsAny(tss_i, tss_j, ignore.strand=FALSE), ]
        tss_i <- sort(reduce(sort(c(tss_i, tss_j[overlapsAny(tss_j, tss_i, ignore.strand=FALSE), ]))))
    }
    tss_i <- sort(reduce(sort(tss_i)))
    tss_all[[i]] <- tss_i
    cat('    ', samples[i], ':\t', length(tss_i), ' TSSs\n', sep = '')
}
tss <- do.call(c, tss_all)
tss <- sort(reduce(sort(tss)))
# What if I extend any TSS smaller than 150, and clip any resulting overlapping TSSs?
tss <- tss[as.character(seqnames(tss)) %in% as.character(1:5), ]
seqlevels(tss) <- as.character(1:5)
tss <- sort(tss)
tss[width(tss) < 150] <- resize(tss[width(tss) < 150], 150, 'center')
tssOvs <- findOverlaps(tss)
tssOvs <- tssOvs[queryHits(tssOvs) != subjectHits(tssOvs)]
if (length(tssOvs)) {
  tssOvs <- tssOvs[seq(1, length(tssOvs), by = 2)]
  if (length(tssOvs)) {
    cat('Adjusting overlapping TSSs ...\n')
    fix_ov_tss <- function(x) {
      wOv <- end(x[1]) - start(x[2])
      mvR <- wOv %/% 2
      mvL <- -(mvR + wOv %% 2)
      x[1] <- shift(x[1], mvL)
      x[2] <- shift(x[2], mvR + 1)
      x
    }
    for (i in seq_len(length(tssOvs))) {
      tss[c(queryHits(tssOvs)[i], subjectHits(tssOvs)[i])] <-
        fix_ov_tss(tss[c(queryHits(tssOvs)[i], subjectHits(tssOvs)[i])])
    }
  }
}
mcols(tss)[['name']] <- paste0('TSS_', 1:length(tss));
cat('Final TSS count: ', length(tss), '\n', sep = '')
export.bed(tss, paste0(\"$OUT_FINAL\", '/tss.final.bed'))
"

echo "Getting raw counts data for merged TSSs ..."

TAGDIRS=(`cat ./tmp_tagdirs`)

if [ ! -f ${OUT_FINAL}/tss.homer.raw.txt ] || [ ! -z $FORCE_QUANT ] ; then
    $HOMER annotatePeaks.pl ${OUT_FINAL}/tss.final.bed \
        tair10 -strand + -fragLength 1 -raw -d ${TAGDIRS[@]} \
        > ${OUT_FINAL}/tss.homer.raw.txt 2> ./tmp_err
fi

echo "Normalizing counts & making bigWigs ..."

$RSCRIPT Rscript -e "
library(edgeR)
library(rtracklayer)

input <- read.table(trimws(\"$INPUT\"), stringsAsFactors=FALSE)

quant <- as.data.frame(readr::read_tsv(trimws(\"${OUT_FINAL}/tss.homer.raw.txt\")))
quant_m <- as.matrix(quant[, 20:ncol(quant)])
colnames(quant_m) <- basename(gsub(' Tag Count .+', '', colnames(quant_m)))
rownames(quant_m) <- quant[[1]]
readr::write_tsv(as.data.frame(cbind(TSS = rownames(quant_m), quant_m)),
    \"${OUT_FINAL}/tss.final.raw.txt\")

y <- DGEList(quant_m, group = gsub('_cs[0-9]$', '', colnames(quant_m)))
y <- calcNormFactors(y, method = 'TMMwsp')

nf <- y[['samples']]
nf[['sample']] <- rownames(nf)
nf[['final.factors']] <- nf[['norm.factors']] * nf[['lib.size']]
nf[['mult.per.million']] <- (nf[['final.factors']] / 1000000)^-1
readr::write_tsv(nf, \"${OUT_FINAL}/norm_factors.txt\")

y_cpm <- cpm(y)
readr::write_tsv(as.data.frame(cbind(TSS = rownames(y_cpm), y_cpm)),
    \"${OUT_FINAL}/tss.final.cpm.txt\")

SL <- c(30427671, 19698289, 23459830, 18585056, 26975502)

samples <- paste0(input[[1]], '_cs', input[[2]])
samples_bg <- paste0(input[[4]], '/', samples)

for (i in seq_along(samples)) {
    cat('    Working on: ', samples[i], '\n', sep = '') 
    bg_p <- suppressWarnings(import(paste0(samples_bg[i], '.raw.pos.bedGraph.gz')))
    bg_n <- suppressWarnings(import(paste0(samples_bg[i], '.raw.neg.bedGraph.gz')))
    bg_p <- bg_p[as.character(seqnames(bg_p)) %in% as.character(1:5), ]
    bg_n <- bg_n[as.character(seqnames(bg_n)) %in% as.character(1:5), ]
    mcols(bg_p)[['score']] <- mcols(bg_p)[['score']] * nf[samples[i], 'mult.per.million']
    mcols(bg_n)[['score']] <- -abs(mcols(bg_n)[['score']]) * nf[samples[i], 'mult.per.million']
    seqlevels(bg_p) <- as.character(1:5)
    seqlevels(bg_n) <- as.character(1:5)
    seqlengths(bg_p) <- SL
    seqlengths(bg_n) <- SL
    export.bw(bg_p, paste0(\"$OUT_BW/\", samples[i], '.rpm.pos.bw'))
    export.bw(bg_n, paste0(\"$OUT_BW/\", samples[i], '.rpm.neg.bw'))
}
"

echo "All done."

rm ./tmp_err

